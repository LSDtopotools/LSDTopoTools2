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
#include <fstream>
#include <math.h>
#include <iostream>
#include "LSDStatsTools.hpp"
#include "LSDParticle.hpp"
using namespace std;

#ifndef LSDParticle_CPP
#define LSDParticle_CPP

// Sorting compiling problems with MSVC
#ifdef _WIN32
#ifndef M_PI
extern double M_PI;
#endif
#endif


const double one_min_exp_neg_2 = 1-exp(-2);
const double one_min_exp_neg_5 = 1-exp(-5);

// default constructor
void LSDParticle::create()
 {
  Type = 0;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
 }

// start with a specific type
void LSDParticle::create( int StartType )
 {
  Type = StartType;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
 }


void LSDParticle::create( int StartType, double StartAge )
 {
  Type = StartType;
  CellIndex = -1;
  Age = StartAge;
 }

void LSDParticle::create( int StartType, int StartCI,  double StartAge, double StartOSLage, double StartxLoc, double StartdLoc)
 {
  Type = StartType;
  CellIndex = StartCI;
  Age = StartAge;
  OSLage = StartOSLage;
  xLoc = StartxLoc;
  dLoc = StartdLoc;
 }

 void LSDParticle::create( int StartType, int StartCI,  double StartAge, 
              double StartOSLage, double StartxLoc, double StartyLoc, double StartdLoc)
 {
  Type = StartType;
  CellIndex = StartCI;
  Age = StartAge;
  OSLage = StartOSLage;
  xLoc = StartxLoc;
  yLoc = StartyLoc;
  dLoc = StartdLoc;
 }

void LSDParticle::create(int StartType, double StartxLoc, double StartdLoc)
 {
  Type = StartType;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
  xLoc = StartxLoc;
  dLoc = StartdLoc;
 }

LSDParticle& LSDParticle::operator=(const LSDParticle& rhs)
 {
  if (&rhs != this)
   {
    create(rhs.getType(),rhs.getCellIndex(),rhs.getAge(),rhs.getOSLage(), rhs.getxLoc(),rhs.getyLoc(),rhs.getdLoc());
   }
  return *this;
 }

std::ostream& operator<<(std::ostream& os, const LSDParticle& tP)
 {
  os <<   tP.getType() << " " << tP.getCellIndex() << " " << tP.getAge() << " "
       << tP.getOSLage() << " " << tP.getxLoc() << " " << tP.getyLoc() 
       << " " << tP.getdLoc();
  return os;
 }

void LSDParticle::incrementAge(double dt)
 {
  if (Age < 0)
   Age = dt;
  else
   Age += dt;
  if (OSLage > 0)
   OSLage += dt;
 }

void LSDParticle::setCellIndex(int tempCI)
{
  CellIndex = tempCI;
}

void LSDParticle::OSLexpose()
 {
   OSLage = 0;
 }

 void LSDParticle::SoilAgeExpose()
{
  Age = 0;
}

// update the x location 
void LSDParticle::update_xLoc(double new_xLoc)
{
  xLoc = new_xLoc;
}

// update the y location
void LSDParticle::update_yLoc(double new_yLoc)
{
  yLoc = new_yLoc;
}


void LSDParticle::displaceReflect(double dx,double dd,double h, double dt)
 {
  xLoc += dx;

  //std::cout << "tPart.cpp LINE 77 dx  is: " << dd << std::endl;
  double dOld = dLoc;
  double dNew = dLoc+dd;
  if (dNew > h)
   {
    //std::cout << "tPart.cpp LINE 84 dNew  is: " << dNew << " and dd is: " << dd
    //          << " and dOld is: "<< dOld << std::endl;
    dLoc = 2*h - dd - dOld;
    //std::cout << "tPart.cpp LINE 86 dLoc  is: " << dLoc << std::endl;
   }
  else if (dNew <= 0)
   {
    dLoc = -(dd+dOld);
    OSLage = (dd+dOld)*dt/dd;
   }
  else
   dLoc = dNew;

  //if (dLoc > h)
  // std::cout << "tPart.cpp LINE 86 dLoc  is: " << dLoc << std::endl;
 }


// this test to see if the particle is still within the x location < lambda
int  LSDParticle::test_domain(double lambda)
 {
  int td;
  if (xLoc >= 0 && xLoc <= lambda)
   td = 1;
  else
   {
    td = -1;
    OSLage = -9999;
    Age = -9999;
    xLoc = -9999;
    dLoc = -9999;
   }
  return td;
 }


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDCRNParticle object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// functions for the CRN loaded tracer particle
void LSDCRNParticle::create()
 {
  double rho_r = 2650;		// in kg/m^3

  Type = 0;
  CellIndex = -1;
  Age = 0;
  OSLage = 0;
  xLoc = 0;			// in metres
  yLoc = 0;     // in metres
  dLoc = 5;			// in metres
  zetaLoc = 100;
  Conc_10Be = 0.0;
  Conc_26Al = 0.0;
  Conc_36Cl = 0.0;
  Conc_14C = 0.0;
  Conc_21Ne = 0.0;
  Conc_3He = 0.0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;
  effective_dLoc = rho_r*dLoc*0.1;	// the 0.1
                    // converts between kg/m^2
                    // to g/cm^2
 }

void LSDCRNParticle::create(int startType, double startxLoc,
              double startzLoc)
{
  double rho_r = 2650;		// in kg/m^3

  Type = startType;
  CellIndex = -1;
  Age = 0;
  OSLage = 0;
  xLoc = startxLoc;			// in meters
  yLoc = 0;
  dLoc = 0;			// in meters
  zetaLoc = startzLoc;
  Conc_10Be = 0.0;
  Conc_26Al = 0.0;
  Conc_36Cl = 0.0;
  Conc_14C = 0.0;
  Conc_21Ne = 0.0;
  Conc_3He = 0.0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;
  effective_dLoc = rho_r*dLoc*0.1;	// the 0.1
                    // converts between kg/m^2
                    // to g/cm^2
 }

void LSDCRNParticle::create(int startType, double startxLoc,
              double startdLoc, double start_effdloc,
              double startzLoc)
{
  Type = startType;
  CellIndex = -1;
  Age = -99;
  OSLage = -99;
  xLoc = startxLoc;			// in meters
  yLoc = 0;
  dLoc = startdLoc;			// in meters
  zetaLoc = startzLoc;
  Conc_10Be = 0.0;
  Conc_26Al = 0.0;
  Conc_36Cl = 0.0;
  Conc_14C = 0.0;
  Conc_21Ne = 0.0;
  Conc_3He = 0.0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;
  effective_dLoc = start_effdloc;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::create(int startType, double startxLoc, double startyLoc,
              double startdLoc, double start_effdloc,
              double startzLoc)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-	            
{
  Type = startType;
  CellIndex = -1;
  Age = -99;
  OSLage = -99;
  xLoc = startxLoc;			// in meters
  yLoc = startyLoc;
  dLoc = startdLoc;			// in meters
  zetaLoc = startzLoc;
  Conc_10Be = 0.0;
  Conc_26Al = 0.0;
  Conc_36Cl = 0.0;
  Conc_14C = 0.0;
  Conc_21Ne = 0.0;
  Conc_3He = 0.0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;
  effective_dLoc = start_effdloc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// a create function for a volume particle
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::create(int startType, int startGSDType,
                            double startxLoc,
                            double startdLoc, double start_effdloc,
                            double startzLoc, double startMass,
                            double startSurfaceArea)
{
  Mass = startMass;					// in kg
  StartingMass = startMass;			// in kg
  SurfaceArea = startSurfaceArea;
                      // in m^2/kg

  Type = startType;
  GSDType = startGSDType;
  CellIndex = -1;
  Age = -99;
  OSLage = -99;
  xLoc = startxLoc;			// in meters
  dLoc = startdLoc;			// in meters
  zetaLoc = startzLoc;
  Conc_10Be = 0.0;
  Conc_26Al = 0.0;
  Conc_36Cl = 0.0;
  Conc_14C = 0.0;
  Conc_21Ne = 0.0;
  Conc_3He = 0.0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;
  effective_dLoc = start_effdloc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-	            
// a create function for a volume particle  with y loc  
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-        
void LSDCRNParticle::create(int startType, int startGSDType, double startxLoc, double startyLoc,
              double startdLoc, double start_effdloc,
              double startzLoc, double startMass,
              double startSurfaceArea)	 
{
  Mass = startMass;                   // in kg
  StartingMass = startMass;           // in kg
  SurfaceArea = startSurfaceArea;
                        // in m^2/kg

  Type = startType;
  GSDType = startGSDType;
  CellIndex = -1;
  Age = -99;
  OSLage = -99;
  xLoc = startxLoc;			// in meters
  yLoc = startyLoc;     // in metres
  dLoc = startdLoc;			// in meters
  zetaLoc = startzLoc;
  Conc_10Be = 0.0;
  Conc_26Al = 0.0;
  Conc_36Cl = 0.0;
  Conc_14C = 0.0;
  Conc_21Ne = 0.0;
  Conc_3He = 0.0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;
  effective_dLoc = start_effdloc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// contains posiution and CRN information
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::create(int startType, double startxLoc,double startzeta_Loc,
                            double startdLoc, double start_effdLoc,
                            double start_C10Be,double start_C14C)
{
  Type = startType;
  CellIndex = -1;
  Age = -99;
  OSLage = -99;
  xLoc = startxLoc;			// in meters
  dLoc = startdLoc;			// in meters
  effective_dLoc = start_effdLoc;
  zetaLoc = startzeta_Loc;
  Conc_10Be = 0.0;
  Conc_26Al = 0.0;
  Conc_36Cl = 0.0;
  Conc_14C = 0.0;
  Conc_21Ne = 0.0;
  Conc_3He = 0.0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This includes some position and CRN information
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::create(int startType, double startxLoc,double startzeta_Loc,
                            double startdLoc, double start_effdLoc,
                            double start_C10Be,double start_C14C, double start_21Ne)
{
  Type = startType;
  CellIndex = -1;
  Age = -99;
  OSLage = -99;
  xLoc = startxLoc;			// in meters
  dLoc = startdLoc;			// in meters
  effective_dLoc = start_effdLoc;
  zetaLoc = startzeta_Loc;
  Conc_10Be = start_C10Be;
  Conc_26Al = 0;
  Conc_36Cl = 0;
  Conc_14C = start_C14C;
  Conc_21Ne = start_21Ne;
  Conc_3He = 0;
  Conc_f7Be = 0.0;
  Conc_f10Be = 0.0;
  Conc_f210Pb = 0.0;
  Conc_f137Cs = 0.0;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This includes all data members for copy constructor
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::create(int startType, int start_GSDType, int startCellIndex, 
        double startAge, double startOSLAge,
        double startxLoc,double startyLoc,double startdLoc, double startefdLoc,
        double startzLoc, double start_C10Be, double start_C26Al,
        double start_C36Cl, double start_C14C,
        double start_C21Ne, double start_C3He,
        double start_Cf7Be, double start_Cf10Be,
        double start_Cf210Pb, double start_Cf137Cs,
        double start_Mass, double start_StartingMass,
        double start_SurfaceArea)
{
  Type = startType;
  CellIndex = startCellIndex;
  Age = startAge;
  OSLage = startOSLAge;
  xLoc = startxLoc;			// in metres
  yLoc = startyLoc;     // in metres
  dLoc = startdLoc;			// in metres
  effective_dLoc = startefdLoc;
  zetaLoc = startzLoc;
  Conc_10Be = start_C10Be;
  Conc_26Al = start_C26Al;
  Conc_36Cl = start_C36Cl;
  Conc_14C = start_C14C;
  Conc_21Ne = start_C21Ne;
  Conc_3He = start_C3He;
  Conc_f7Be = start_Cf7Be;
  Conc_f10Be = start_Cf10Be;
  Conc_f210Pb = start_Cf210Pb;
  Conc_f137Cs = start_Cf137Cs;
  Mass = start_Mass;
  StartingMass = start_StartingMass;
  SurfaceArea =  start_SurfaceArea;
  GSDType = start_GSDType;	
}


LSDCRNParticle& LSDCRNParticle::operator=(const LSDCRNParticle& rhs)
{
  if (&rhs != this)
    {
      create(rhs.getType(), rhs.getGSDType(), rhs.getCellIndex(), 
               rhs.getAge(),rhs.getOSLage(),
               rhs.getxLoc(),rhs.getyLoc(),rhs.getdLoc(),
               rhs.geteffective_dLoc(),rhs.get_zetaLoc(),rhs.getConc_10Be(),
               rhs.getConc_26Al(), rhs.getConc_36Cl(), rhs.getConc_14C(),
               rhs.getConc_21Ne(), rhs.getConc_3He(),
               rhs.getConc_f7Be(), rhs.getConc_f10Be(),
               rhs.getConc_f210Pb(), rhs.getConc_f137Cs(),
               rhs.getMass(), rhs.getStartingMass(),
               rhs.getSurfaceArea());
    }
    return *this;
}

LSDCRNParticle& LSDCRNParticle::operator=(LSDCRNParticle& rhs)
{
  if (&rhs != this)
    {
      create(rhs.getType(), rhs.getGSDType(), rhs.getCellIndex(), 
               rhs.getAge(),rhs.getOSLage(),
               rhs.getxLoc(),rhs.getyLoc(),rhs.getdLoc(),
               rhs.geteffective_dLoc(),rhs.get_zetaLoc(),rhs.getConc_10Be(),
               rhs.getConc_26Al(), rhs.getConc_36Cl(), rhs.getConc_14C(),
               rhs.getConc_21Ne(), rhs.getConc_3He(),
               rhs.getConc_f7Be(), rhs.getConc_f10Be(),
               rhs.getConc_f210Pb(), rhs.getConc_f137Cs(),
               rhs.getMass(), rhs.getStartingMass(),
               rhs.getSurfaceArea());
    }
    return *this;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Unit conversion utilities
// length units are in metres or g/cm^2
// density is in kg/m^3
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::convert_m_to_gpercm2(double l_in_m, double rho)
{
  if(rho < 500)
  {
    cout << "Density should be in kg/m^3, your density is: " << rho << endl;
    cout << "Are you sure you have the correct units?" << endl;
  }
  
  double l_in_gpercm2 = l_in_m*rho*0.1;
  return l_in_gpercm2;
}

double LSDCRNParticle::convert_gpercm2_to_m(double l_in_gpercm2, double rho)
{
  if(rho < 500)
  {
    cout << "Density should be in kg/m^3, your density is: " << rho << endl;
    cout << "Are you sure you have the correct units?" << endl;
  }
  
  double l_in_m = l_in_gpercm2*10/rho;
  return l_in_m;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function just sets the concentration of cosmogenic in situ nuclides
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_cosmo_conc_const(double C_10Be, double C_26Al, double C_36Cl,
                          double C_14C, double C_21Ne, double C_3He)
{
  Conc_10Be = C_10Be;
  Conc_26Al = C_26Al;
  Conc_36Cl = C_36Cl;
  Conc_14C  = C_14C;
  Conc_21Ne = C_21Ne;
  Conc_3He  = C_3He;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function updates the concentration of CRNs in a particle
// the model assumes that during the timestep the change in the
// 'depth' of the particle occurs ofver a constant rate.
// The depth in this case is an equvalent depth...it is linearly
// proportional to depth if overlying density is constant, but
// really represetnt the mass per area above a point in the subsurface
// thus the erosion rate represent a mass removal rate.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_10Be_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  // fist berillium
  double Be_exp = exp(-dt*CRNp.lambda_10Be);
  //double Be_depth_exp = exp(-d_0/Gamma_Be1);

  //cout << "LINE 236 LSDParticle.cpp " << endl;
  //cout << lambda_10Be << " " << Be_exp << endl;
  //cout << "starting conc: " << Conc_10Be << endl;
  //cout << "starting depth: " << effective_dLoc << endl;
  //cout << "erosion_rate: "<< erosion_rate << endl;
  double sum_term = 0;
  for (int i = 0; i<4; i++)
  {
    sum_term+= (CRNp.F_10Be[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
           (exp(dt*erosion_rate/CRNp.Gamma[i])-
            exp(-dt*CRNp.lambda_10Be))/
           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
  //cout << "F_10Be["<<i<<"]: " << CRNp.F_10Be[i] << " and sum_term: " << sum_term <<endl;
  }

  //cout << "and sum term is: " << sum_term << endl;
  
  // Note: Scaling is either done by modifying S_t or by scaling the F factors
  // for the basinwide cosmogenic tool S_t is set to 1 and the scaling is done
  // with the F factors.   
  Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*CRNp.P0_10Be*sum_term;
  //cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}

// this function updates the 10Be concentration if there is a linear
// increase (or decrease) in erosion rate.
void LSDCRNParticle::update_10Be_conc_linear_increase(double dt,double erosion_rate, double alpha, LSDCRNParameters& CRNp)
{
  // fist berillium
  double Be_exp = exp(-dt*CRNp.lambda_10Be);
  //double Be_depth_exp = exp(-d_0/Gamma_Be1);

    //cout << "LINE 236 LSDParticle.cpp " << endl;
    //cout << lambda_10Be << " " << Be_exp << endl;
    //cout << "starting conc: " << Conc_10Be << endl;
    //cout << "starting depth: " << effective_dLoc << endl;
    //cout << "erosion_rate: "<< erosion_rate << endl;
  double sum_term = 0;
  double L_j, M_j, A_j, B_j, D_j;
  double erfi_arg1, erfi_arg2;
  for (int j = 0; j<4; j++)
  {
    A_j = sqrt(0.5*CRNp.Gamma[j]/alpha);
    B_j = sqrt(2*alpha*CRNp.Gamma[j]);
    D_j = sqrt(0.5*alpha/CRNp.Gamma[j]);
    L_j = exp(-( (erosion_rate+CRNp.Gamma[j]*CRNp.lambda_10Be)*(erosion_rate+CRNp.Gamma[j]*CRNp.lambda_10Be)+
            2*alpha*(effective_dLoc+dt*CRNp.Gamma[j]*CRNp.lambda_10Be) )/(B_j*B_j) );
    erfi_arg1 = erosion_rate/B_j + A_j*CRNp.lambda_10Be;
    erfi_arg2 = D_j*dt + erfi_arg1;
    M_j = erfi(erfi_arg2) - erfi(erfi_arg1);

    sum_term+= (CRNp.F_10Be[j]*A_j*L_j*M_j);
  }

    //cout << "and sum term is: " << sum_term << endl;

  Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*CRNp.P0_10Be*sum_term*sqrt(M_PI);
  //cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}

void LSDCRNParticle::update_10Be_conc_neutron_only(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{

  double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
  double Be_exp = exp(-dt*CRNp.lambda_10Be);

  double sum_term = 0;
  sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
           (exp(dt*erosion_rate/Gamma_neutron)-
            exp(-dt*CRNp.lambda_10Be))/
           (erosion_rate+Gamma_neutron*CRNp.lambda_10Be);

  Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*Be_exp*CRNp.P0_10Be*sum_term;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This updates the 10Be concentration if erosion is steady using the full
// 4 attenuation depth model of Vermeesch (2007)
// The erosion rate should be in g/cm^2/yr
//
// NOTE This produces far less muon production than CRONUS
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_10Be_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
  double this_term;
  double sum_term1 = 0;
  //double sum_term2 = 0;
  
  double spall_tot = 0; 
  double muon_tot = 0;
  
  for (int i = 0; i<4; i++)
  {
    //cout << "ThickSF: " <<  CRNp.F_10Be[i]*(1-exp(-effective_dLoc/CRNp.Gamma[i]))*
    //                        (CRNp.Gamma[i]/effective_dLoc) << endl;
    //cout << "A: " <<  (CRNp.lambda_10Be+erosion_rate/CRNp.Gamma[i]) << endl;
  
    //sum_term1+= (CRNp.F_10Be[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
    //           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
    //sum_term1+= (CRNp.F_10Be[i]*(1-exp(-effective_dLoc/CRNp.Gamma[i]))*
    //            (CRNp.Gamma[i]/effective_dLoc))/
    //            (CRNp.lambda_10Be+erosion_rate/CRNp.Gamma[i]);
    //cout << "Thick/A:" << (CRNp.F_10Be[i]*(1-exp(-effective_dLoc/CRNp.Gamma[i]))*
    //            (CRNp.Gamma[i]/effective_dLoc))/
    //            (CRNp.lambda_10Be+erosion_rate/CRNp.Gamma[i]) << endl;
    this_term = (exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.F_10Be[i]*CRNp.Gamma[i])/
                 (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
    sum_term1 += this_term;
    if(i == 0)
    {
      spall_tot+=this_term;
    }             
    else
    {
      muon_tot+=this_term;
    }
  }

  //cout << "effective_dLoc: " << effective_dLoc  << endl;
  //cout << "erosion rate: " << erosion_rate << endl;
  

  //cout << "and sum term is: " << sum_term1 << " " << sum_term2 << endl;
  //double Pref =  CRNp.S_t*CRNp.P0_10Be;
  //cout << "Pref is: " << Pref << endl;
  //cout << "Scaling is: " << CRNp.S_t << " P0 is: " << CRNp.P0_10Be 
  //     << " and Pref is: " << CRNp.S_t*CRNp.P0_10Be << endl;
  Conc_10Be = CRNp.S_t*CRNp.P0_10Be*sum_term1;
  spall_tot = CRNp.S_t*CRNp.P0_10Be*spall_tot;
  muon_tot =  CRNp.S_t*CRNp.P0_10Be*muon_tot;
  //cout << "Line 717, Conc 10Be is: " << Conc_10Be << " from spallation: " << spall_tot
  //     << " and muons: " << muon_tot << endl;
  
  //cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// A depth integrated version of the steady state equation. 
// This is really only appropriate for basinwide calculations because 
// it doens't really represent the CRN concentration of an individual particle
// The arguments also override the effective depth of the particle
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_10Be_SSfull_depth_integrated(double erosion_rate, 
                                  LSDCRNParameters& CRNp,
                                  double top_eff_depth, double bottom_eff_depth)
{

  //cout << "LINE 748, doing full depth integration" << endl;

  // first, check if the inputs are working properly
  if (top_eff_depth > bottom_eff_depth)
  {
    cout << "LSDParticle line 753, your effective depths for integration are reversed" << endl;
    cout << "Reversing the two depths. Check your inputs!" << endl;
    double temp_eff_depth = bottom_eff_depth;
    bottom_eff_depth = top_eff_depth;
    top_eff_depth = temp_eff_depth;
  }

  // If the top effective depth and the bottom effective depth are the same, then run the standard
  // concentration for particles eroding from the top of this pixel
  if (top_eff_depth == bottom_eff_depth)
  {
    effective_dLoc = top_eff_depth;
    update_10Be_SSfull(erosion_rate, CRNp);
  }
  else  // If they arent the same, this calculates the depth average concentrations of particles between the top and bottom effective depths.  
  {
    double this_term;
    double sum_term1 = 0;
  
    double spall_tot = 0; 
    double muon_tot = 0;
    
    for (int i = 0; i<4; i++)
    { 
      // get the individual terms
      this_term = ( ( exp(-top_eff_depth/CRNp.Gamma[i])
                     -exp(-bottom_eff_depth/CRNp.Gamma[i]) )*
                    CRNp.F_10Be[i]*CRNp.Gamma[i]*CRNp.Gamma[i])/
                   (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
      sum_term1 += this_term;
      
      // add the terms to track muon and spalling production
      if(i == 0)
      {
        spall_tot+=this_term;
      }             
      else
      {
        muon_tot+=this_term;
      }
      
    }
    
    // get the concentration
    Conc_10Be = (1/(bottom_eff_depth-top_eff_depth))*CRNp.S_t*CRNp.P0_10Be*sum_term1;
    spall_tot = (1/(bottom_eff_depth-top_eff_depth))*CRNp.S_t*CRNp.P0_10Be*spall_tot;
    muon_tot =  (1/(bottom_eff_depth-top_eff_depth))*CRNp.S_t*CRNp.P0_10Be*muon_tot;

    
    
  }

  // print some stuff to screen for bug checking  
  vector<bool> cosmo_flags(4,false);
  cosmo_flags[0] = true;

  //cout << "\n\n\nLSDParticle line 797, getting cosmo for erosion rate: " << erosion_rate <<endl;
  //cout << "Top depth: " << top_eff_depth << " bottom depth: " << bottom_eff_depth << endl;
  //cout << "The parameter values are: " << endl;
  //CRNp.print_parameters_to_screen(cosmo_flags);
  //cout << "Conc 10Be: " << Conc_10Be << endl << endl << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function uses newton iteration to get an apparent erosion rate
// for full muons. 
// It is intended to mimic the cosmocalc erosion rate finder
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParticle::apparent_erosion_10Be_COSMOCALC(double rho, LSDCRNParameters& CRNp,
                                    double scaling_factor, string Muon_scaling,
                                    double top_eff_depth, double bottom_eff_depth)
{
  double Target_conc = Conc_10Be;
  
  // get the initial guess
  effective_dLoc = top_eff_depth;
  double erate_guess = apparent_erosion_10Be_neutron_only(rho, CRNp);
  double eff_erate_guess = convert_m_to_gpercm2(erate_guess,rho);
  
  //cout << "Concentration: " << Conc_10Be << " rho: " << rho << " neutron erate " 
  //     << erate_guess << "m/yr, effective: " << eff_erate_guess << " g/cm^2/yr" << endl;

  // set the scaling
  // reset scaling parameters. This is necessary since the F values are
  // reset for local scaling
  if (Muon_scaling == "Schaller" )
  {
    CRNp.set_Schaller_parameters();
  }
  else if (Muon_scaling == "Braucher" )
  {
    CRNp.set_Braucher_parameters();
  }
  else if (Muon_scaling == "Granger" )
  {
    CRNp.set_Granger_parameters();
  }
  else if (Muon_scaling == "newCRONUS" )
  {
    CRNp.set_newCRONUS_parameters();
  }
  else
  {
    cout << "You didn't set the muon scaling." << endl
         << "Options are Schaller, Braucher, newCRONUS, and Granger." << endl
         << "You chose: " << Muon_scaling << endl
         << "Defaulting to Braucher et al (2009) scaling" << endl;
    CRNp.set_Braucher_parameters();     
  }
  
  // scale the F values
  vector<bool> nuclide_scaling_switches(4,false);
  nuclide_scaling_switches[0] = true;       // this scales for 10Be
  CRNp.scale_F_values(scaling_factor,nuclide_scaling_switches);
  
  // now do Newton Iteration to find the correct erosion rate
  // convert to  g/cm^2/yr
  eff_erate_guess = 0.1*erate_guess*rho;
  
  // now using this as the initial guess, use Newton-Raphson to zero in on the
  // correct erosion rate
  double eff_e_new = eff_erate_guess; // the erosion rate upon which we iterate
  double eff_e_change;                // the change in erosion rate between iterations
  double tolerance = 1e-10;           // tolerance for a change in the erosion rate
                                      // between Newton-Raphson iterations
  double eff_e_displace = 1e-6;       // A small displacment in the erosion rate used
                                      // to calculate the derivative
  double N_this_step;                 // the concentration of the nuclide reported this step
  double N_displace;                  // the concentration at the displaced erosion rate
  double N_derivative;                // dN/de derivative for Newton-Raphson
  double f_x;                         // the function being tested by newton raphson
  double f_x_displace;                // the displaced function (for calculating the derivative)
  
  do
  {
    update_10Be_SSfull_depth_integrated(eff_e_new, CRNp,
                                           top_eff_depth,  bottom_eff_depth);
    N_this_step = Conc_10Be;
    
    update_10Be_SSfull_depth_integrated(eff_e_new+eff_e_displace, CRNp,
                                           top_eff_depth,  bottom_eff_depth); 
    N_displace = Conc_10Be;
    
    
    // print to screen:
    //cout << "eff_e_new = " << eff_e_new << " in cm/kyr: " << eff_e_new*1e6/rho 
    //     << " 10Be: " << N_this_step << " atoms/yr" << endl;
    
    f_x =  N_this_step-Target_conc;
    f_x_displace =  N_displace-Target_conc;
    
    N_derivative = (f_x_displace-f_x)/eff_e_displace;
      
    if(N_derivative != 0)
    {
      
      eff_e_new = eff_e_new-f_x/N_derivative;
      
      // check to see if the difference in erosion rates meet a tolerance
      eff_e_change = f_x/N_derivative;
      //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
    }
    else
    {
      eff_e_change = 0;
    }
  
  } while(fabs(eff_e_change) > tolerance);
  
  vector<double> erosion_rate_vec;
  erosion_rate_vec.push_back(eff_e_new);
  erosion_rate_vec.push_back(eff_e_new*10/rho);
  
  Conc_10Be = Target_conc;
  return erosion_rate_vec;

}                                    
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_10Be_neutron_only(double rho, LSDCRNParameters& CRNp)
{
  // a few constants, all computed from Vermeesh 2007
  double Gamma_neutron=CRNp.Gamma[0];     // in g/cm^2
  double app_eff_eros;                    // in g/cm2/yr
  double app_eros;                        // in m/yr
  
  double exp_term = exp(-effective_dLoc/Gamma_neutron);

  //cout << "LSDCRNParticle line 741, S_t: " << CRNp.neutron_S_t << " P0: " << CRNp.P0_10Be 
  //     << " total spallation: " << CRNp.neutron_S_t*CRNp.P0_10Be << endl 
  //     << " and the F_0: " << CRNp.F_10Be[0] << endl;

  app_eff_eros = Gamma_neutron*(exp_term*CRNp.neutron_S_t*CRNp.P0_10Be
                                /Conc_10Be-CRNp.lambda_10Be);
  app_eros = app_eff_eros*10/rho;                              
  return app_eros;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


void LSDCRNParticle::update_26Al_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  double Al_exp = exp(-dt*CRNp.lambda_26Al);

  double sum_term = 0;
  for (int i = 0; i<4; i++)
  {
    sum_term+= (CRNp.F_26Al[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
           (exp(dt*erosion_rate/CRNp.Gamma[i])-
            exp(-dt*CRNp.lambda_26Al))/
           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
  }

  Conc_26Al = Conc_26Al*Al_exp +  CRNp.S_t*Al_exp*CRNp.P0_26Al*sum_term;
}



void LSDCRNParticle::update_26Al_conc_neutron_only(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2
  double Al_exp = exp(-dt*CRNp.lambda_26Al);

  double sum_term = 0;
  sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
           (exp(dt*erosion_rate/Gamma_neutron)-
            exp(-dt*CRNp.lambda_26Al))/
           (erosion_rate+Gamma_neutron*CRNp.lambda_26Al);

  Conc_26Al = Conc_26Al*Al_exp +  CRNp.S_t*Al_exp*CRNp.P0_26Al*sum_term;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_26Al_neutron_only(double rho, LSDCRNParameters& CRNp)
{
  // a few constants, all computed from Vermeesh 2007
  double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
  double app_eff_eros;                    // in g/cm2/yr
  double app_eros;                        // in m/yr
  
  double exp_term = exp(-effective_dLoc/Gamma_neutron);

  app_eff_eros = Gamma_neutron*(exp_term*CRNp.neutron_S_t*CRNp.P0_26Al
                                /Conc_26Al-CRNp.lambda_26Al);
  app_eros = app_eff_eros*10/rho;                              
  return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function uses newton iteration to get an apparent erosion rate
// for full muons. 
// It is intended to mimic the cosmocalc erosion rate finder
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParticle::apparent_erosion_26Al_COSMOCALC(double rho, LSDCRNParameters& CRNp,
                                    double scaling_factor, string Muon_scaling,
                                    double top_eff_depth, double bottom_eff_depth)
{
  double Target_conc = Conc_26Al;
  
  // get the initial guess
  effective_dLoc = top_eff_depth;
  double erate_guess = apparent_erosion_26Al_neutron_only(rho, CRNp);
  double eff_erate_guess = convert_m_to_gpercm2(erate_guess,rho);
  
  //cout << "Concentration: " << Conc_26Al << " rho: " << rho << " neutron erate " 
  //     << erate_guess << "m/yr, effective: " << eff_erate_guess << " g/cm^2/yr" << endl;

  // set the scaling
  // reset scaling parameters. This is necessary since the F values are
  // reset for local scaling
  if (Muon_scaling == "Schaller" )
  {
    CRNp.set_Schaller_parameters();
  }
  else if (Muon_scaling == "Braucher" )
  {
    CRNp.set_Braucher_parameters();
  }
  else if (Muon_scaling == "Granger" )
  {
    CRNp.set_Granger_parameters();
  }
  else if (Muon_scaling == "newCRONUS" )
  {
    CRNp.set_newCRONUS_parameters();
  }
  else
  {
    cout << "You didn't set the muon scaling." << endl
         << "Options are Schaller, Braucher, newCRONUS, and Granger." << endl
         << "You chose: " << Muon_scaling << endl
         << "Defaulting to Braucher et al (2009) scaling" << endl;
    CRNp.set_Braucher_parameters();     
  }
  
  // scale the F values
  vector<bool> nuclide_scaling_switches(4,false);
  nuclide_scaling_switches[1] = true;       // this scales for 26Al
  CRNp.scale_F_values(scaling_factor,nuclide_scaling_switches);
  
  // now do Newton Iteration to find the correct erosion rate
  // convert to  g/cm^2/yr
  eff_erate_guess = 0.1*erate_guess*rho;
  
  // now using this as the initial guess, use Newton-Raphson to zero in on the
  // correct erosion rate
  double eff_e_new = eff_erate_guess; // the erosion rate upon which we iterate
  double eff_e_change;                // the change in erosion rate between iterations
  double tolerance = 1e-10;           // tolerance for a change in the erosion rate
                                      // between Newton-Raphson iterations
  double eff_e_displace = 1e-6;       // A small displacment in the erosion rate used
                                      // to calculate the derivative
  double N_this_step;                 // the concentration of the nuclide reported this step
  double N_displace;                  // the concentration at the displaced erosion rate
  double N_derivative;                // dN/de derivative for Newton-Raphson
  double f_x;                         // the function being tested by newton raphson
  double f_x_displace;                // the displaced function (for calculating the derivative)
  
  do
  {
    update_26Al_SSfull_depth_integrated(eff_e_new, CRNp,
                                           top_eff_depth,  bottom_eff_depth);
    N_this_step = Conc_26Al;
    
    update_26Al_SSfull_depth_integrated(eff_e_new+eff_e_displace, CRNp,
                                           top_eff_depth,  bottom_eff_depth); 
    N_displace = Conc_26Al;
    
    
    // print to screen:
    //cout << "eff_e_new = " << eff_e_new << " in cm/kyr: " << eff_e_new*1e6/rho 
    //     << " 26Al: " << N_this_step << " atoms/yr" << endl;
    
    f_x =  N_this_step-Target_conc;
    f_x_displace =  N_displace-Target_conc;
    
    N_derivative = (f_x_displace-f_x)/eff_e_displace;
      
    if(N_derivative != 0)
    {
      
      eff_e_new = eff_e_new-f_x/N_derivative;
      
      // check to see if the difference in erosion rates meet a tolerance
      eff_e_change = f_x/N_derivative;
      //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
    }
    else
    {
      eff_e_change = 0;
    }
  
  } while(fabs(eff_e_change) > tolerance);
  
  vector<double> erosion_rate_vec;
  erosion_rate_vec.push_back(eff_e_new);
  erosion_rate_vec.push_back(eff_e_new*10/rho);
  
  Conc_26Al = Target_conc;
  return erosion_rate_vec;

}                                    
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Gets steady concentration
// The erosion rate should be in g/cm^2/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_26Al_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
  double this_sum;
  double sum_term = 0;
  double spall_tot = 0; 
  double muon_tot = 0;
  for (int i = 0; i<4; i++)
  {
    //sum_term+= (CRNp.F_26Al[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
    //     (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
    this_sum = (exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.F_26Al[i]*CRNp.Gamma[i])/
                 (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
    sum_term+=this_sum;
    if(i == 0)
    {
      spall_tot+=this_sum;
      //cout << "Stot: " << spall_tot << endl;
    }             
    else
    {
      muon_tot+=this_sum;
      //cout << "Mtot: " << muon_tot << endl;
    }                 
  }

  Conc_26Al = CRNp.S_t*CRNp.P0_26Al*sum_term;
  spall_tot = CRNp.S_t*CRNp.P0_26Al*spall_tot;
  muon_tot =  CRNp.S_t*CRNp.P0_26Al*muon_tot;
  cout << "Line 717, Conc 26Al is: " << Conc_26Al << " from spallation: " << spall_tot
       << " and muons: " << muon_tot << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// A depth integrated version of the steady state equation. 
// This is really only appropriate for basinwide calculations because 
// it doens't really represent the CRN concentration of an individual particle
// The arguments also override the effective depth of the particle
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_26Al_SSfull_depth_integrated(double erosion_rate, 
                                  LSDCRNParameters& CRNp,
                                  double top_eff_depth, double bottom_eff_depth)
{

  // first, check if the inputs are working properly
  if (top_eff_depth > bottom_eff_depth)
  {
    cout << "LSDParticle line 753, your effective depths for integration are reversed" << endl;
    cout << "Reversing the two depths. Check your inputs!" << endl;
    double temp_eff_depth = bottom_eff_depth;
    bottom_eff_depth = top_eff_depth;
    top_eff_depth = temp_eff_depth;
  }

  if (top_eff_depth == bottom_eff_depth)
  {
    effective_dLoc = top_eff_depth;
    update_26Al_SSfull(erosion_rate, CRNp);
  }
  else
  {
    double this_term;
    double sum_term1 = 0;
  
    double spall_tot = 0; 
    double muon_tot = 0;
    
    for (int i = 0; i<4; i++)
    { 
      // get the individual terms
      this_term = ( ( exp(-top_eff_depth/CRNp.Gamma[i])
                     -exp(-bottom_eff_depth/CRNp.Gamma[i]) )*
                    CRNp.F_26Al[i]*CRNp.Gamma[i]*CRNp.Gamma[i])/
                   (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
      sum_term1 += this_term;
      
      // add the terms to track muon and spalling production
      if(i == 0)
      {
        spall_tot+=this_term;
      }             
      else
      {
        muon_tot+=this_term;
      }
      
    }
    
    // get the concentration
    Conc_26Al = (1/(bottom_eff_depth-top_eff_depth))*CRNp.S_t*CRNp.P0_26Al*sum_term1;
    spall_tot = (1/(bottom_eff_depth-top_eff_depth))*CRNp.S_t*CRNp.P0_26Al*spall_tot;
    muon_tot =  (1/(bottom_eff_depth-top_eff_depth))*CRNp.S_t*CRNp.P0_26Al*muon_tot;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_14C_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  double C_exp = exp(-dt*CRNp.lambda_14C);

  double sum_term = 0;
  for (int i = 0; i<4; i++)
  {
    sum_term+= (CRNp.F_14C[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
           (exp(dt*erosion_rate/CRNp.Gamma[i])-
            exp(-dt*CRNp.lambda_14C))/
           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_14C);
  }

  Conc_14C = Conc_14C*C_exp +  CRNp.S_t*C_exp*CRNp.P0_14C*sum_term;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_14C_conc_neutron_only(double dt,double erosion_rate, 
                                                   LSDCRNParameters& CRNp)
{
  double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

  double C_exp = exp(-dt*CRNp.lambda_14C);

  double sum_term = 0;
  sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
           (exp(dt*erosion_rate/Gamma_neutron)-
            exp(-dt*CRNp.lambda_14C))/
           (erosion_rate+Gamma_neutron*CRNp.lambda_14C);

  Conc_14C = Conc_14C*C_exp +  CRNp.S_t*C_exp*CRNp.P0_14C*sum_term;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_14C_neutron_only(double rho, LSDCRNParameters& CRNp)
{
  // a few constants, all computed from Vermeesh 2007
  double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
  double app_eff_eros;                    // in g/cm2/yr
  double app_eros;                        // in m/yr
  
  double exp_term = exp(-effective_dLoc/Gamma_neutron);

  app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_14C
                                /Conc_14C-CRNp.lambda_14C);
  app_eros = app_eff_eros*10/rho;                              
  return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


// The erosion rate should be in g/cm^2/yr
void LSDCRNParticle::update_14C_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
  double sum_term = 0;
  for (int i = 0; i<4; i++)
  {
    sum_term+= (CRNp.F_14C[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
               (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_14C);
  }

  Conc_14C = CRNp.S_t*CRNp.P0_14C*sum_term;
}

void LSDCRNParticle::update_36Cl_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  double Cl_exp = exp(-dt*CRNp.lambda_36Cl);

  double sum_term = 0;
  for (int i = 0; i<4; i++)
  {
    sum_term+= (CRNp.F_36Cl[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
           (exp(dt*erosion_rate/CRNp.Gamma[i])-
            exp(-dt*CRNp.lambda_36Cl))/
           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_36Cl);
  }

  Conc_36Cl = Conc_36Cl*Cl_exp +  CRNp.S_t*Cl_exp*CRNp.P0_36Cl*sum_term;
}

void LSDCRNParticle::update_36Cl_conc_neutron_only(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

  double Cl_exp = exp(-dt*CRNp.lambda_36Cl);

  double sum_term = 0;
  sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
           (exp(dt*erosion_rate/Gamma_neutron)-
            exp(-dt*CRNp.lambda_36Cl))/
           (erosion_rate+Gamma_neutron*CRNp.lambda_36Cl);


  Conc_36Cl = Conc_36Cl*Cl_exp +  CRNp.S_t*Cl_exp*CRNp.P0_36Cl*sum_term;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_36Cl_neutron_only(double rho, LSDCRNParameters& CRNp)
{
  // a few constants, all computed from Vermeesh 2007
  double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
  double app_eff_eros;                    // in g/cm2/yr
  double app_eros;                        // in m/yr
  
  double exp_term = exp(-effective_dLoc/Gamma_neutron);

  app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_36Cl
                                /Conc_36Cl-CRNp.lambda_36Cl);
  app_eros = app_eff_eros*10/rho;                              
  return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



// The erosion rate should be in g/cm^2/yr
void LSDCRNParticle::update_36Cl_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
  double sum_term = 0;
  for (int i = 0; i<4; i++)
  {
    sum_term+= (CRNp.F_36Cl[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
               (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_36Cl);
  }

  Conc_36Cl = CRNp.S_t*CRNp.P0_36Cl*sum_term;
}

void LSDCRNParticle::update_21Ne_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  double Gamma_neutron= CRNp.Gamma[0];					// in g/cm^2
  if (erosion_rate == 0)
  {
    Conc_21Ne = Conc_21Ne +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*CRNp.P0_21Ne*dt;
  }
  else
  {
    Conc_21Ne = Conc_21Ne +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_21Ne*
               (exp(dt*erosion_rate/Gamma_neutron) - 1)/erosion_rate;
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_21Ne(double rho, LSDCRNParameters& CRNp)
{
  // a few constants, all computed from Vermeesh 2007
  double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
  double app_eff_eros;                    // in g/cm2/yr
  double app_eros;                        // in m/yr
  
  double exp_term = exp(-effective_dLoc/Gamma_neutron);

  app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_21Ne
                                /Conc_21Ne);
  app_eros = app_eff_eros*10/rho;                              
  return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



// The erosion rate should be in g/cm^2/yr
void LSDCRNParticle::update_21Ne_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
  double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

  //cout << endl << endl <<"BEFORE Conc_21Ne: " << Conc_21Ne << endl;
  if (erosion_rate*erosion_rate < 0.0000000001)
  {
    Conc_21Ne = 0;
  }
  else
  {
    Conc_21Ne = CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_21Ne/erosion_rate;
  }
  //cout << "AFTER Conc_21Ne: " << Conc_21Ne << endl;
}

void LSDCRNParticle::update_3He_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
  double Gamma_neutron = CRNp.Gamma[0];	// in g/cm^2
  Conc_3He = Conc_3He +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_3He*
               (exp(dt*erosion_rate/Gamma_neutron) - 1)/erosion_rate;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_3He(double rho, LSDCRNParameters& CRNp)
{
  // a few constants, all computed from Vermeesh 2007
  double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
  double app_eff_eros;                    // in g/cm2/yr
  double app_eros;                        // in m/yr
  
  double exp_term = exp(-effective_dLoc/Gamma_neutron);

  app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_3He
                                /Conc_3He);
  app_eros = app_eff_eros*10/rho;                              
  return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


// The erosion rate should be in g/cm^2/yr
void LSDCRNParticle::update_3He_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
  double Gamma_neutron = CRNp.Gamma[0];	// in g/cm^2

  if (erosion_rate*erosion_rate < 0.0000000001)
  {
    Conc_3He = 0;
  }
  else
  {
    Conc_3He = Conc_3He +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_3He/erosion_rate;
 }
}

void LSDCRNParticle::update_all_CRN(double dt, double erosion_rate, LSDCRNParameters& CRNp)
{
  //cout << "LINE 445 LSDCRNParticle.cpp updating CRN" << endl;
  update_10Be_conc(dt,erosion_rate, CRNp);
  update_26Al_conc(dt,erosion_rate, CRNp);
  update_36Cl_conc(dt,erosion_rate, CRNp);
  update_14C_conc(dt,erosion_rate, CRNp);
  update_21Ne_conc(dt,erosion_rate, CRNp);
  update_3He_conc(dt,erosion_rate, CRNp);
}


// The erosion rate should be in g/cm^2/yr
void LSDCRNParticle::update_all_CRN_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
  //cout << "LINE 445 LSDCRNParticle.cpp updating CRN" << endl;
  update_10Be_SSfull(erosion_rate, CRNp);
  update_26Al_SSfull(erosion_rate, CRNp);
  update_36Cl_SSfull(erosion_rate, CRNp);
  update_14C_SSfull(erosion_rate, CRNp);
  update_21Ne_SSfull(erosion_rate, CRNp);
  update_3He_SSfull(erosion_rate, CRNp);
}

void LSDCRNParticle::update_all_CRN_neutron_only(double dt, double erosion_rate, LSDCRNParameters& CRNp)
{
  update_10Be_conc_neutron_only(dt,erosion_rate, CRNp);
  update_26Al_conc_neutron_only(dt,erosion_rate, CRNp);
  update_36Cl_conc_neutron_only(dt,erosion_rate, CRNp);
  update_14C_conc_neutron_only(dt,erosion_rate, CRNp);
  update_21Ne_conc(dt,erosion_rate, CRNp);
  update_3He_conc(dt,erosion_rate, CRNp);
}

// this updates fallout radionuclides
// NOTE!!!!
// the untis of M_supply are in atoms/cm^2
// and conc is in atoms per gram
// we convert rho into g/cm^3
// adn deltad into cm
// withing the LSDCRNParticle object
// the units of k_f10Be here are cm^2/g
void LSDCRNParticle::update_fallout10Be_simple_density(double dt, double M_supply_surface,
          double rho_skg, double k_f10Be, double deltad_m, LSDCRNParameters& CRNp)
{
  // first find which depth interval the particle is in
  int depth_interval = int(dLoc/deltad_m);
  double d_top = double(depth_interval)*deltad_m*100;
  double d_bottom = double(depth_interval+1)*deltad_m*100;
  double deltad = deltad_m*100;
  // the factor of 100 is to convert to cm

  // convert density to g/cm^3
  double rho_s = rho_skg/1000;

  // get the cutoff depth
  double cutoff_depth = 5/(rho_s*k_f10Be);

  if (dLoc*100 > cutoff_depth)
  {
    Conc_f10Be += - Conc_f10Be*CRNp.lambda_10Be;
  }
  else
  {
    Conc_f10Be += dt*M_supply_surface*( exp(-k_f10Be*rho_s*d_top) -exp(-k_f10Be*rho_s*d_bottom) )/
                  (deltad*rho_s*one_min_exp_neg_5) - Conc_f10Be*CRNp.lambda_10Be;
  }
}

// this updates fallout radionuclides
// NOTE!!!!
// the units of M_supply are in atoms/cm^2
// and conc is in atoms per gram
// we convert rho into g/cm^3
// adn deltad into cm
// within the LSDCRNParticle object
// the units of k_f10Be here are cm^2/g
// chi_f10Be is the fraction of shallow fallout 10Be
void LSDCRNParticle::update_fallout10Be_simple_density_2exp(double dt, double M_supply_surface,
                              double rho_skg, double k1_f10Be, double k2_f10Be, 
                              double chi_f10Be, double deltad_m, LSDCRNParameters& CRNp)
{
  // first find which depth interval the particle is in
  int depth_interval = int(dLoc/deltad_m);
  double d_top = double(depth_interval)*deltad_m*100;
  double d_bottom = double(depth_interval+1)*deltad_m*100;
  double deltad = deltad_m*100;
  // the factor of 100 is to convert to cm

  // convert density to g/cm^3
  double rho_s = rho_skg/1000;

  // get the cutoff depth for k1
  double cutoff_depth1 = 5/(rho_s*k1_f10Be);
  double cutoff_depth2 = 5/(rho_s*k2_f10Be);

  if (dLoc*100 > cutoff_depth2)
  {
    Conc_f10Be += - Conc_f10Be*CRNp.lambda_10Be;
  }
  if (dLoc*100 < cutoff_depth2 && dLoc*100 > cutoff_depth1)
  {
    Conc_f10Be += dt*M_supply_surface*(1-chi_f10Be)*( exp(-k2_f10Be*rho_s*d_top) -exp(-k2_f10Be*rho_s*d_bottom) )/
                  (deltad*rho_s*one_min_exp_neg_5) - Conc_f10Be*CRNp.lambda_10Be;
  }
  else
  {
    Conc_f10Be += dt*M_supply_surface*(1-chi_f10Be)*( exp(-k2_f10Be*rho_s*d_top) -exp(-k2_f10Be*rho_s*d_bottom) )/
                   (deltad*rho_s*one_min_exp_neg_5) +
                  dt*M_supply_surface*chi_f10Be*( exp(-k1_f10Be*rho_s*d_top) -exp(-k1_f10Be*rho_s*d_bottom) )/
                   (deltad*rho_s*one_min_exp_neg_5) -
                  Conc_f10Be*CRNp.lambda_10Be;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function resets the depth and effective depth
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_depths(double d, double ed)
{
  dLoc=d;
  effective_dLoc=ed;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function resets the depth and calculates the effective depth
// based on a 1 layer model
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::calculate_effective_depth_one_layer(double d, double rho)
{
  dLoc=d;
  
  // now get the effective depth
  // the 0.1 is because effective depth is in g.cm^2
  // depth is in metres, rho is in kg/m^3
  // so a factor of 0.1 is to convert kg/m^2 to g/cm^2
  effective_dLoc = rho*0.1*d;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function changes the effective depth (i.e., the shielding depth)
// for a constant rate of erosion
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::erode_mass_only(double dt, double mass_erosion_rate)
{
  effective_dLoc-= mass_erosion_rate*dt;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function changes the effective depth (i.e., the shielding depth)
// for erosion that is changing linearly in time
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::erode_mass_only_linear_increase(double dt, double mass_erosion_rate, double alpha)
{
  effective_dLoc-= 0.5*dt*(alpha*dt+2*mass_erosion_rate);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function simply resets zeta
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_zetaLoc(double new_zeta)
{
  zetaLoc = new_zeta;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function resets the zeta locations using an updated
// surface elevation that preserves the d and effective d locations
// it is only for use with a model that has no soil
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_zetaLoc_with_new_surface(double new_zeta)
{
  zetaLoc = new_zeta-dLoc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// !!!FUNCTIONS FOR SCALINF
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// for Deselits 2006 scaling
// -------- function definition for denominator of DZ2006 Eqns. 4 and 7 -------
//that is, the integral of the expression for Beta
// Modified from Greg Balco's code:
//  http://hess.ess.washington.edu/math
//
// x integrating variable
// Rc is cutoff rigidity (GV)
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::intOfBeta(double x, double Rc)
{

  // coefficients from 2006 paper
  double n = 1.0177e-2;
  double alpha = 1.0207e-1;
  double k = -3.9527e-1;
  double a0 = 8.5236e-6;
  double a1 = -6.3670e-7;
  double a2 = -7.0814e-9;
  double a3 = -9.9182e-9;
  double a4 = 9.9250e-10;
  double a5 = 2.4925e-11;
  double a6 = 3.8615e-12;
  double a7 = -4.8194e-13;
  double a8 = -1.5371e-14;

  double out = n*x/(1 + exp(-alpha*(pow(Rc,-k)))) + (0.5)*(a0+a1*Rc+a2*Rc*Rc)*(x*x)+
    (1/3.0)*(a3+a4*Rc+a5*Rc*Rc)*(x*x*x) + (0.25)*(a6+a7*Rc+a8*Rc*Rc)*(x*x*x*x);

  return out;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The deselits 2006 sclaing
// Modified from Greg Balco's code:
//  http://hess.ess.washington.edu/math
//
// h is atmospheric pressure (hPa)
// Rc is cutoff rigidity (GV)
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::desilets2006sp(double h,double Rc)
{
  double x = h*1.019716;

  // enforce rigidity knee at 2 GV;
  if (Rc < 2.0)
  {
    Rc = 2.0;
  }

  // reference atmospheric depth at SL
  double sld = 1033.0;

  //get altitude scaling factor
  // first get attenuation length
  double L = (sld - x)/(intOfBeta(sld,Rc) - intOfBeta(x,Rc));

  // apply attenuation length
  double fofx = exp((sld-x)/L);

  // get latitude scaling factor
  double alpha = 10.275;
  double k = 0.9615;

  double fofRc = 1 - exp(-alpha*( pow(Rc,-k)));

  // now the total scaling factor
  double out = fofRc*fofx;
  return out;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets the Dunai Scaling
// Modified from Greg Balco's code:
//  http://hess.ess.washington.edu/math
//
// h is atmospheric pressure (hPa)
// Rc is cutoff rigidity (GV)
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double LSDCRNParticle::dunai2001sp(double h, double Rc)
{
  double delz = (1013.25-h)*1.019716;

  // constants
  double A = 0.5221;
  double B = -1.7211;
  double C = 0.3345;
  double X = 4.2822;
  double Y = 0.4952;

  double a = 17.183;
  double b = 2.060;
  double c = 5.9164;
  double x = 2.2964;
  double y = 130.11;

  // sea level scaling factor
  double N1030 = Y + A/( pow( (1 + exp(-(Rc-X)/B)),C));

  // attenuation length
  double L = y + a/( pow((1 + exp(-(Rc-x)/b)),c));

  //total scaling factor
  double out = N1030*exp(delz/L);

  return out;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets the Lifton 2006 Scaling
// Modified from Greg Balco's code:
//  http://hess.ess.washington.edu/math
//
// h is atmospheric pressure (hPa)
// Rc is cutoff rigidity (GV)
// S solar modulation factor (nondimensional, see source paper)
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double LSDCRNParticle::lifton2006sp(double h,double Rc,double S)
{
  // convert pressure to atmospheric depth
  double X = h*1.019716;

  // flatten low rigidities. The value of 1.907 comes from the source paper.
  if (Rc < 1.907)
  {
    Rc = 1.907;
  }

  //define constants
  vector<double> c;
  c.push_back(1.8399);
  c.push_back(-1.1854e2);
  c.push_back(-4.9420e-2);
  c.push_back(8.0139e-1);
  c.push_back(1.2708e-4);
  c.push_back(9.4647e-1);
  c.push_back(-3.2208e-2);
  c.push_back(1.2688);

  double t1 = c[0]*log(X*S);
  double t2 = -S*exp( (c[1]*S)/(   pow((Rc + 5.0*S),(2.0*S)) ) );
  double t3 = c[2]*(pow(X,c[3]));
  double t4 = c[4]*( pow( (Rc + 4.0*S)*X,c[5]));
  double t5 = c[6]*( pow( (Rc + 4.0*S),c[7]));

  double out = exp(t1 + t2 + t3 + t4 + t5);
  return out;
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
double LSDCRNParticle::stone2000sp(double lat,double P, double Fsp)
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
// This implements a cutoff-rigidity rather than latitude based scaling
// scheme based on the Lal spallation polynomials. For use in
// paleomagnetically-corrected exposure-age calculations. 
//
//   h = atmospheric pressure (hPa)
//   Rc = cutoff rigidity (GV)
//
// Original version Written by Greg Balco -- UW Cosmogenic Nuclide Lab
// balcs@u.washington.edu
// March, 2007
// Part of the CRONUS-Earth online calculators: 
//      http://hess.ess.washington.edu/math
//
// IMPORTANT: This (and the stonesp version) is probably the best scaling method!
// See https://cosmognosis.wordpress.com/2014/01/07/high-altitude-low-latitude-calibration-sites-i/
//
// Updated for c++ by Simon Mudd
// 05/12/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::stone2000Rcsp(double h, double Rc)
{

  if (Rc > 21)
  {
    cout << "Your cutoff rigidity is greater than 21GV. " << endl;
    cout << "Defaulting to 21 GV" << endl;
    Rc = 21;
  }

  // Build the scaling factor = f(rigidity) function up to 14.9 GV
  vector<double> ilats_degree(7,0.0);
  vector<double> ilats_radians(7,0.0);
  ilats_degree[0] = 0;
  ilats_radians[0] = ilats_degree[0]*M_PI/180.0;
  ilats_degree[1] = 10;
  ilats_radians[1] = ilats_degree[1]*M_PI/180.0;
  ilats_degree[2] = 20;
  ilats_radians[2] = ilats_degree[2]*M_PI/180.0;
  ilats_degree[3] = 30;
  ilats_radians[3] = ilats_degree[3]*M_PI/180.0;
  ilats_degree[4] = 40;
  ilats_radians[4] = ilats_degree[4]*M_PI/180.0;
  ilats_degree[5] = 50;
  ilats_radians[5] = ilats_degree[5]*M_PI/180.0;
  ilats_degree[6] = 60;
  ilats_radians[6] = ilats_degree[6]*M_PI/180.0;            

  // Convert latitude to rigidity using Elsasser formula (from Sandstrom)
  // Rigidity Rc = Rc(0)cos^4(latitude)
  // where Rc(0) = rigidity at equator, that is, 14.9 GV
  vector<double> iRcs(8,0.0);
  for(int i = 0; i<7; i++)
  {
    
    iRcs[i] = 14.9*(cos(ilats_radians[i])*cos(ilats_radians[i])*cos(ilats_radians[i])*cos(ilats_radians[i]));  
    cout << "Ilats rad["<<i+1<<"]: " <<ilats_radians[i] << " and iRcs: " << iRcs[i] <<  endl;
  }
 
  // Now Spallogenic production at index rigidities;
  // Constants from Table 1 of Stone(2000)
  vector<double> a(7,0.0);
  vector<double> b(7,0.0);
  vector<double> c(7,0.0);
  vector<double> d(7,0.0);
  vector<double> e(7,0.0);
  
  a[0] = 31.8518;
  a[1] = 34.3699;
  a[2] = 40.3153;
  a[3] = 42.0983;
  a[4] = 56.7733;
  a[5] = 69.0720;
  a[6] = 71.8733;
  
  b[0] = 250.3193;
  b[1] = 258.4759;
  b[2] = 308.9894;
  b[3] = 512.6857;
  b[4] = 649.1343;
  b[5] = 832.4566;
  b[6] = 863.1927;
  
  c[0] = -0.083393;
  c[1] = -0.089807;
  c[2] = -0.106248;
  c[3] = -0.120551;
  c[4] = -0.160859;
  c[5] = -0.199252;
  c[6] = -0.207069;
  
  d[0] = 7.4260e-5;
  d[1] = 7.9457e-5;
  d[2] = 9.4508e-5;
  d[3] = 1.1752e-4;
  d[4] = 1.5463e-4;
  d[5] = 1.9391e-4;
  d[6] = 2.0127e-4;
  
  e[0] = -2.2397e-8;
  e[1] = -2.3697e-8; 
  e[2] = -2.8234e-8;
  e[3] = -3.8809e-8;
  e[4] = -5.0330e-8;
  e[5] = -6.3653e-8;
  e[6] = -6.6043e-8;

  // Apply Eqn. (2) of Stone (2000);
  vector<double> sf(8,0.0);
  sf[0] = a[0] + ( b[0]*exp(h/(-150.0)) + c[0]*h + d[0]*h*h + e[0]*h*h*h);
  sf[1] = a[1] + ( b[1]*exp(h/(-150.0)) + c[1]*h + d[1]*h*h + e[1]*h*h*h);
  sf[2] = a[2] + ( b[2]*exp(h/(-150.0)) + c[2]*h + d[2]*h*h + e[2]*h*h*h);
  sf[3] = a[3] + ( b[3]*exp(h/(-150.0)) + c[3]*h + d[3]*h*h + e[3]*h*h*h);
  sf[4] = a[4] + ( b[4]*exp(h/(-150.0)) + c[4]*h + d[4]*h*h + e[4]*h*h*h);
  sf[5] = a[5] + ( b[5]*exp(h/(-150.0)) + c[5]*h + d[5]*h*h + e[5]*h*h*h);
  sf[6] = a[6] + ( b[6]*exp(h/(-150.0)) + c[6]*h + d[6]*h*h + e[6]*h*h*h);
  
  // Extend to zero rigidity - scaling factor does not change from that at 60
  // degrees --
  iRcs[7] = 0; 
  sf[7] = sf[6];

  // Extend to 21 GV by fitting a log-log line to the latitude 0-20 values,
  // i.e. where rigidity is greater than 10 GV. According to Quemby and Wenk, 
  // as summarized  in Sandstrom, log(rigidity) vs. log(nucleon intensity) 
  // ought to be linear above 10 GV. Note that this is speculative, but 
  // relatively unimportant, as the approximation is pretty much only used 
  // for low latitudes for a short time in the Holocene. 

  // convert some vect to floats so the linear regression works
  // this is a bit stupid but I haven't bothered to template things
  vector<float> ln_iRcs_float(3,0.0);
  vector<float> ln_sf_float(3,0.0);
  vector<float> resid(3,0.0);
  for(int i = 0; i<3; i++)
  {
    ln_iRcs_float[i] = float(log(iRcs[i]));
    ln_sf_float[i] = float(log(sf[i]));
  }
  
  vector<float> linfit = simple_linear_regression(ln_iRcs_float,ln_sf_float,resid);
  
  vector<double> new_sf(14,0.0);
  vector<double> new_iRcs(14,0.0);
  
  cout << "Linfit m: " << linfit[0] << endl;
  cout << "Linfit b: " << linfit[1] << endl;
  
  for(int i = 0; i< 6; i++)
  {
    new_sf[i] = exp( log(sf[0] + double(linfit[0])*(log(21.0-double(i)) - log(iRcs[0]))));
    new_iRcs[i] = 21.0-double(i);
  }
  for(int i = 0; i<8; i++)
  {
    new_sf[i+6] = sf[i];
    new_iRcs[i+6] = iRcs[i];
  }
  
  
  for(int i = 0; i<14; i++)
  {
    cout << "i: " << i << " iRcs: " << new_iRcs[i] << " sf: " << new_sf[i] <<  endl;
  }

  cout << "Entering interp1D" << endl;

  // Interpolate, return
  double out = interp1D_unordered(new_iRcs,new_sf,Rc);
  cout << "Scaling is: " << out << endl;
  
  return out;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function computes the length scaling
// It is from the CRONUS calculator but includes a flag that uses the
// Lambda data members
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::thickness_scaling_factor(LSDCRNParameters& LSDCRNP, bool use_CRONUS)
{
  double this_Gamma;
  //double gamma_metric;
  double thick_scale;
  if(effective_dLoc > 0)
  {
    this_Gamma = LSDCRNP.get_spallation_attenuation_length(use_CRONUS);
    
    //cout << "gamma: " << this_Gamma << " and eff d: "<< effective_dLoc << endl;
    
    thick_scale = (this_Gamma/effective_dLoc)*
                        (1- exp(-effective_dLoc/this_Gamma));
    //cout << "Thick scale is: " << thick_scale << endl;                    
  }
  else
  {
    thick_scale = 1;
  }
  return thick_scale;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is the initial guess of erosion rate
// that repicates the initial guesss component of the CRONUS 
// erosion rate function
// The erosion returned is in g/cm^2/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParticle::CRONUS_initial_guess(LSDCRNParameters& LSDCRNP, double pressure,
                                 double lat, double N_10Be, double N_26Al, 
                                 double topo_scale, double snow_scale)
{

  // get the other shielding correction
  double other_shielding = topo_scale*snow_scale;
  
  // get the muon production
  double test_elev = 0.0;
  vector<double> muon_prod 
             = LSDCRNP.calculate_muon_production_CRONUS(test_elev, pressure);
  vector<double> Prefs_st = LSDCRNP.get_Stone_Pref();
  
  // get the Stone scaling
  double Fsp = 1.0;     // for the initial guess we don't adjust Fsp (as in CRONUS)
  double stoneP = stone2000sp(lat, pressure, Fsp);
  //cout << "LINE 1718 Lat: " << lat << " pressure: " << pressure << " stone: " << stoneP << endl;
  
  // retrieve the pre-scaling factors
  double P_ref_St_10 = Prefs_st[0];
  double P_ref_St_26 = Prefs_st[1];
  

  // get the scaling from thickness
  bool use_CRONUS = true;
  double this_thickSF = thickness_scaling_factor(LSDCRNP, use_CRONUS);
  
  double P_temp_10 = (P_ref_St_10*stoneP*this_thickSF*other_shielding)
                       + muon_prod[0] + muon_prod[2];
  double P_temp_26 = (P_ref_St_26*stoneP*this_thickSF*other_shielding)
                       + muon_prod[1] + muon_prod[3];                       

  //cout << "10Be fast: "  << muon_prod[0] << " 10Be neg" << muon_prod[2] << endl;
  //cout << "P temp 10: " << P_temp_10 << endl;


  double Gamma = LSDCRNP.get_spallation_attenuation_length(use_CRONUS);
  vector<double> decay_coeff = LSDCRNP.get_decay_coefficients(use_CRONUS);
  
  //cout << "Gamma: " << Gamma << " N_10Be" << N_10Be << " lambda: " << decay_coeff[0] << endl;
  
  
  double E_lal_10 = 0.0;
  double E_lal_26 = 0.0;
  if(N_10Be > 0)
  {
    E_lal_10 = Gamma*(P_temp_10/N_10Be - decay_coeff[0]);
  }
  if(N_26Al > 0)
  {
    E_lal_26 = Gamma*(P_temp_26/N_26Al - decay_coeff[1]);
  }

  vector<double> erosion_guess(2,0.0);
  erosion_guess[0] =  E_lal_10;
  erosion_guess[1] =  E_lal_26;
  
  return erosion_guess;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// CRONUS Get erosion
// This function emulates the CroNUS calculator
// get_al_be_erosion.m
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParticle::CRONUS_get_Al_Be_erosion(LSDCRNParameters& LSDCRNP, double pressure,
                      double lat, double rho, double N_10Be, double N_26Al, 
                      double sample_del10, double sample_del26,
                      double topo_scale, double snow_scale)
{
  // some erosion parameters that will be determined by this function 
  // erate_consts holds all the information that replicates the CRONUS calculator
  // e.g., erosion rate, production rates and uncertainties. 
  double eff_e_10Be = 0;
  double eff_e_26Al = 0;
  vector<double> erate_consts(12,0.0);
  
  bool use_CRONUS = true;
  
  // first scale the thickness
  double thickSF = thickness_scaling_factor(LSDCRNP, use_CRONUS);
  cout << "thickSF is: " << thickSF << endl;
  
  // Now get the initial guess
  vector<double> initial_guess = CRONUS_initial_guess(LSDCRNP, pressure, lat, 
                                         N_10Be, N_26Al, topo_scale, snow_scale);
  
  // variable for holding the number of atoms per gram
  double N10_this_step,N10_displace;
  //double N10;
  double N26_this_step,N26_displace;
  //double N26;
  double N_derivative;        // derivative of the function describing N as a fxn of erosion
  double e_displace = 0.0001;   // how far one moves erosion rate to calculate derivative
  double e_new,e_change;      // the erosion rate in the Newton-Raphson iteration
  double tolerance = 1e-10;
  double f_x;                 // we call this the 'function' where we need to find
                              // a root. It is just the N atoms minus the target atoms
  double f_x_displace;                           
  
  // get some parameters for Stone production
  vector<double> Prefs_st = LSDCRNP.get_Stone_Pref();
  
  // get the Stone scaling
  double Fsp = 1.0;     // for the initial guess we don't adjust Fsp (as in CRONUS)
  double stoneP = stone2000sp(lat, pressure, Fsp);
  //cout << "LINE 1718 Lat: " << lat << " pressure: " << pressure << " stone: " << stoneP << endl;
  
  // retrieve the pre-scaling factors
  double P_ref_St_10 = Prefs_st[0];
  double P_ref_St_26 = Prefs_st[1];
  //cout << "P ref stone 10: " <<  P_ref_St_10 <<  " P_ref_St_26: " << P_ref_St_26 << endl;
  
  // precalculate the P_mu vectors
  vector<double> z_mu;
  vector<double> P_mu_z_10Be;
  vector<double> P_mu_z_26Al;
  LSDCRNP.get_CRONUS_P_mu_vectors(pressure, effective_dLoc, z_mu, 
                                  P_mu_z_10Be,P_mu_z_26Al);
                                  
  // get the scaled production rates
  double P_sp_10Be = P_ref_St_10*stoneP*topo_scale*snow_scale;
  double P_sp_26Al = P_ref_St_26*stoneP*topo_scale*snow_scale;
  
  //cout << "P_sp_10Be: " << P_sp_10Be << " P_sp_26Al: " << P_sp_26Al << endl;
  
  if(N_10Be > 0)
  {
    // now enter search for the correct erosion rate. 
    // first do 10Be
    // get the initial concentrations
    e_new = initial_guess[0];
    //cout << endl << endl << "10Be initial guess is: " << e_new << endl;
    // now iterate using newton-raphson
    do
    {
      
      // get the new values
      CRONUS_calculate_N_forward(e_new, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_this_step, N26_this_step);
       f_x =  N10_this_step-N_10Be;                  
                          
      // now get the derivative
      CRONUS_calculate_N_forward(e_new+e_displace, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_displace, N26_displace);
      f_x_displace =  N10_displace-N_10Be;                    
                          
      N_derivative = (f_x_displace-f_x)/e_displace;
      
      //cout << "N10: " << N10_this_step << " + displace: " << N10_displace  << endl;
      //cout << "fx: " <<  f_x << " f_x_displace: " << f_x_displace << " fprimex: " << N_derivative << endl;
      
      if(N_derivative != 0)
      {
        e_new = e_new-f_x/N_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        e_change = f_x/N_derivative;
        //cout << "Change is: " << e_change << " and erosion rate is: " << e_new << endl;
      }
      else
      {
        e_change = 0;
      }
    
    } while(fabs(e_change) > tolerance);
    eff_e_10Be = e_new; 
  }

  // now enter search for the correct erosion rate. 
  // now do 26Al
  
  if(N_26Al > 0)
  {
    // get the initial concentrations
    e_new = initial_guess[1];
    //cout << endl << endl << "26Al initial guess is: " << e_new << endl;
    // now iterate using newton-raphson
    do
    {
      
      // get the new values
      CRONUS_calculate_N_forward(e_new, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_this_step, N26_this_step);
       f_x =  N26_this_step-N_26Al;                  
                          
      // now get the derivative
      CRONUS_calculate_N_forward(e_new+e_displace, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_displace, N26_displace);
      f_x_displace =  N26_displace-N_26Al;                    
                          
      N_derivative = (f_x_displace-f_x)/e_displace;
      
      //cout << "N10: " << N10_this_step << " + displace: " << N10_displace  << endl;
      //cout << "fx: " <<  f_x << " f_x_displace: " << f_x_displace << " fprimex: " << N_derivative << endl;
      
      if(N_derivative != 0)
      {
        e_new = e_new-f_x/N_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        e_change = f_x/N_derivative;
        //cout << "Change is: " << e_change << " and erosion rate is: " << e_new << endl;
      }
      else
      {
        e_change = 0;
      }
    
    } while(fabs(e_change) > tolerance);
    eff_e_26Al = e_new;
  }
  
  vector<double> uncertainties = 
           CRONUS_error_propagation(pressure, LSDCRNP, thickSF,rho,N_10Be, N_26Al,
                      sample_del10, sample_del26, z_mu,  P_mu_z_10Be, P_mu_z_26Al,
                      P_sp_10Be, P_sp_26Al, eff_e_10Be, eff_e_26Al);
  
  cout << "10Be atoms from spallation: " << uncertainties[8] 
       << " and muons: " << uncertainties[9]  << endl;
  cout << "26Al atoms from spallation: " << uncertainties[10] 
       << " and muons: " << uncertainties[11]  << endl;       
  
  
  erate_consts[0] = eff_e_10Be;
  erate_consts[1] = eff_e_10Be*1.0e7/rho;
  erate_consts[2] = uncertainties[1];
  erate_consts[3] = uncertainties[0];
  erate_consts[4] = uncertainties[3];
  erate_consts[5] = P_sp_10Be;
  erate_consts[6] = eff_e_26Al;
  erate_consts[7] = eff_e_26Al*1.0e7/rho;
  erate_consts[8] = uncertainties[5];
  erate_consts[9] = uncertainties[4];
  erate_consts[10] = uncertainties[7];
  erate_consts[11] = P_sp_26Al; 
  
  return erate_consts;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// CRONUS Get erosion
// This function emulates the CroNUS calculator
// get_al_be_erosion.m
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParticle::CRONUS_get_Al_Be_erosion_modified_production(LSDCRNParameters& LSDCRNP, double pressure,
                      double lat, double rho, double N_10Be, double N_26Al, 
                      double sample_del10, double sample_del26,
                      double topo_scale, double snow_scale,
                      double Spallation_frac, double P_mu_frac)
{
  // some erosion parameters that will be determined by this function 
  // erate_consts holds all the information that replicates the CRONUS calculator
  // e.g., erosion rate, production rates and uncertainties. 
  double eff_e_10Be = 0;
  double eff_e_26Al = 0;
  vector<double> erate_consts(12,0.0);
  
  bool use_CRONUS = true;
  
  // first scale the thickness
  double thickSF = thickness_scaling_factor(LSDCRNP, use_CRONUS);
  cout << "thickSF is: " << thickSF << endl;
  
  // Now get the initial guess
  vector<double> initial_guess = CRONUS_initial_guess(LSDCRNP, pressure, lat, 
                                         N_10Be, N_26Al, topo_scale, snow_scale);
  
  // variable for holding the number of atoms per gram
  double N10_this_step,N10_displace;
  //double N10;
  double N26_this_step,N26_displace;
  //double N26;
  double N_derivative;        // derivative of the function describing N as a fxn of erosion
  double e_displace = 0.0001;   // how far one moves erosion rate to calculate derivative
  double e_new,e_change;      // the erosion rate in the Newton-Raphson iteration
  double tolerance = 1e-10;
  double f_x;                 // we call this the 'function' where we need to find
                              // a root. It is just the N atoms minus the target atoms
  double f_x_displace;                           
  
  // get some parameters for Stone production
  vector<double> Prefs_st = LSDCRNP.get_Stone_Pref();
  
  // get the Stone scaling
  double Fsp = 1.0;     // for the initial guess we don't adjust Fsp (as in CRONUS)
  double stoneP = stone2000sp(lat, pressure, Fsp);
  //cout << "LINE 1718 Lat: " << lat << " pressure: " << pressure << " stone: " << stoneP << endl;
  
  // retrieve the pre-scaling factors
  double P_ref_St_10 = Prefs_st[0]*Spallation_frac;
  double P_ref_St_26 = Prefs_st[1]*Spallation_frac;
  //cout << "P ref stone 10: " <<  P_ref_St_10 <<  " P_ref_St_26: " << P_ref_St_26 << endl;
  
  // precalculate the P_mu vectors
  vector<double> z_mu;
  vector<double> P_mu_z_10Be;
  vector<double> P_mu_z_26Al;
  LSDCRNP.get_CRONUS_P_mu_vectors(pressure, effective_dLoc, z_mu, 
                                  P_mu_z_10Be,P_mu_z_26Al);
  
  int n_pmu = (P_mu_z_10Be.size());
  for (int i = 0; i<n_pmu; i++)
  {
    P_mu_z_10Be[i] = P_mu_z_10Be[i]*P_mu_frac;
    P_mu_z_26Al[i] = P_mu_z_26Al[i]*P_mu_frac;
  }
  
  // get the scaled production rates
  double P_sp_10Be = P_ref_St_10*stoneP*topo_scale*snow_scale;
  double P_sp_26Al = P_ref_St_26*stoneP*topo_scale*snow_scale;
  
  //cout << "P_sp_10Be: " << P_sp_10Be << " P_sp_26Al: " << P_sp_26Al << endl;
  
  if(N_10Be > 0)
  {
    // now enter search for the correct erosion rate. 
    // first do 10Be
    // get the initial concentrations
    e_new = initial_guess[0];
    //cout << endl << endl << "10Be initial guess is: " << e_new << endl;
    // now iterate using newton-raphson
    do
    {
      
      // get the new values
      CRONUS_calculate_N_forward(e_new, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_this_step, N26_this_step);
       f_x =  N10_this_step-N_10Be;                  
                          
      // now get the derivative
      CRONUS_calculate_N_forward(e_new+e_displace, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_displace, N26_displace);
      f_x_displace =  N10_displace-N_10Be;                    
                          
      N_derivative = (f_x_displace-f_x)/e_displace;
      
      //cout << "N10: " << N10_this_step << " + displace: " << N10_displace  << endl;
      //cout << "fx: " <<  f_x << " f_x_displace: " << f_x_displace << " fprimex: " << N_derivative << endl;
      
      if(N_derivative != 0)
      {
        e_new = e_new-f_x/N_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        e_change = f_x/N_derivative;
        //cout << "Change is: " << e_change << " and erosion rate is: " << e_new << endl;
      }
      else
      {
        e_change = 0;
      }
    
    } while(fabs(e_change) > tolerance);
    eff_e_10Be = e_new; 
  }

  // now enter search for the correct erosion rate. 
  // now do 26Al
  
  if(N_26Al > 0)
  {
    // get the initial concentrations
    e_new = initial_guess[1];
    //cout << endl << endl << "26Al initial guess is: " << e_new << endl;
    // now iterate using newton-raphson
    do
    {
      
      // get the new values
      CRONUS_calculate_N_forward(e_new, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_this_step, N26_this_step);
       f_x =  N26_this_step-N_26Al;                  
                          
      // now get the derivative
      CRONUS_calculate_N_forward(e_new+e_displace, LSDCRNP,z_mu, 
                          P_mu_z_10Be, P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                          N10_displace, N26_displace);
      f_x_displace =  N26_displace-N_26Al;                    
                          
      N_derivative = (f_x_displace-f_x)/e_displace;
      
      //cout << "N10: " << N10_this_step << " + displace: " << N10_displace  << endl;
      //cout << "fx: " <<  f_x << " f_x_displace: " << f_x_displace << " fprimex: " << N_derivative << endl;
      
      if(N_derivative != 0)
      {
        e_new = e_new-f_x/N_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        e_change = f_x/N_derivative;
        //cout << "Change is: " << e_change << " and erosion rate is: " << e_new << endl;
      }
      else
      {
        e_change = 0;
      }
    
    } while(fabs(e_change) > tolerance);
    eff_e_26Al = e_new;
  }
  
  vector<double> uncertainties = 
           CRONUS_error_propagation(pressure, LSDCRNP, thickSF,rho,N_10Be, N_26Al,
                      sample_del10, sample_del26, z_mu,  P_mu_z_10Be, P_mu_z_26Al,
                      P_sp_10Be, P_sp_26Al, eff_e_10Be, eff_e_26Al);
  
  cout << "10Be atoms from spallation: " << uncertainties[8] 
       << " and muons: " << uncertainties[9]  << endl;
  cout << "26Al atoms from spallation: " << uncertainties[10] 
       << " and muons: " << uncertainties[11]  << endl;       
  
  
  erate_consts[0] = eff_e_10Be;
  erate_consts[1] = eff_e_10Be*1.0e7/rho;
  erate_consts[2] = uncertainties[1];
  erate_consts[3] = uncertainties[0];
  erate_consts[4] = uncertainties[3];
  erate_consts[5] = P_sp_10Be;
  erate_consts[6] = eff_e_26Al;
  erate_consts[7] = eff_e_26Al*1.0e7/rho;
  erate_consts[8] = uncertainties[5];
  erate_consts[9] = uncertainties[4];
  erate_consts[10] = uncertainties[7];
  erate_consts[11] = P_sp_26Al; 
  
  return erate_consts;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function produces a screen report similar to that in CRONUS calculator
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::CRONUS_screen_report(vector<double> erate_consts)
{

  cout << endl << endl << endl;
  cout << "||==============================================================================================||" << endl;
  cout << "|| Reporting 10Be results for non-time dependent Stone scaling.                                 ||" << endl;
  cout << "||______________________________________________________________________________________________||" << endl;
  cout << "|| Muon Prod \t|| Internal\t||Erate    \t||Erate    \t|| External\t|| Spallation\t||" << endl;
  cout << "|| rate      \t|| Uncert. \t||g/cm^2/yr\t||m/Myr    \t|| Uncert. \t|| prod rate \t||" << endl;
  cout << "|| atoms/g/yr\t|| m/Myr   \t||         \t||         \t|| m/Myr   \t|| atoms/g/yr \t||" << endl;
  cout << "||-----------\t||---------\t||---------\t||---------\t||---------\t||------------\t||" << endl;
  cout << "|| " << erate_consts[4] << "\t|| " << erate_consts[2] << "\t|| " << erate_consts[0] 
       << " || " << erate_consts[1] << "\t|| " << erate_consts[3] << "\t|| " 
       << erate_consts[5] << "\t||" << endl;
  cout << "||==============================================================================================||" << endl;
  cout << "||==============================================================================================||" << endl;  
  cout << "|| Reporting 26Al results for non-time dependent Stone scaling.                                 ||" << endl;
  cout << "||______________________________________________________________________________________________||" << endl;
  cout << "|| Muon Prod \t|| Internal\t||Erate    \t||Erate    \t|| External\t|| Spallation \t||" << endl;
  cout << "|| rate      \t|| Uncert. \t||g/cm^2/yr\t||m/Myr    \t|| Uncert. \t|| prod rate  \t||" << endl;
  cout << "|| atoms/g/yr\t|| m/Myr   \t||         \t||         \t|| m/Myr   \t|| atoms/g/yr \t||" << endl;
  cout << "||-----------\t||---------\t||---------\t||---------\t||---------\t||------------\t||" << endl;
  cout << "|| " << erate_consts[10] << "\t|| " << erate_consts[8] << "\t|| " 
       << erate_consts[6] << "\t|| " << erate_consts[7] << "\t|| " 
       << erate_consts[9] << "\t|| " << erate_consts[11] << "\t||" <<endl;
  cout << "||==============================================================================================||" << endl;  
  cout << endl << endl;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function returns the number of atoms given an effective erosion rate
// (this is the erosion rate in g/cm^2/yr)
// It is used in an optimisation loop (equivalent to fzero in matlab)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::CRONUS_calculate_N_forward(double effective_erosion_rate, 
                            LSDCRNParameters& LSDCRNP,
                            vector<double>& z_mu, vector<double>& P_mu_z_10Be, 
                            vector<double>& P_mu_z_26Al, double thickSF, 
                            double P_sp_10Be, double P_sp_26Al,
                            double& N_Be10, double& N_Al26)
{                            
  // get the initial fluxes from muons and spallation
  double Be10_mu_N;
  double Al26_mu_N;
  LSDCRNP.integrate_muon_flux_for_erosion(effective_erosion_rate,z_mu, P_mu_z_10Be,
                           P_mu_z_26Al, Be10_mu_N, Al26_mu_N);
                           
  //cout << "Initial muon production guess 10Be: " << Be10_mu_N << " Al26: " << Al26_mu_N << endl;

  double Be10_sp_N;
  double Al26_sp_N;
  LSDCRNP.integrate_nonTD_spallation_flux_for_erosion(effective_erosion_rate,thickSF,
                           P_sp_10Be, P_sp_26Al,Be10_sp_N, Al26_sp_N);

  N_Be10 = Be10_sp_N+Be10_mu_N;
  N_Al26 = Al26_sp_N+Al26_mu_N;
  //cout << "Initial spallation production guess 10Be: " << Be10_sp_N << " Al26: " << Al26_sp_N << endl;                                                                                                                             
  
  //============================================================================
  // test these against the nuclides from analytical muons (e.g., based
  // on the granger or Schaller equations)
  //double erate_m_yr = convert_gpercm2_to_m(initial_guess[0]);
  //LSDCRNP.set_Neutron_only_parameters();
  //LSDCRNP.set_CRONUS_stone_parameters();
  //LSDCRNP.set_scaling(stoneP, topo_scale, snow_scale);
  //update_10Be_SSfull(initial_guess[0],LSDCRNP);
  //update_26Al_SSfull(initial_guess[0],LSDCRNP);
  //cout << "Analytical solution 10: " << Conc_10Be << " and 26: " << Conc_26Al << endl;
  //============================================================================
  // I'VE DONE THIS AND CRONUS CALC PREDICTS WAY MORE MUON PRODUCTION (sometimes
  // as much as 50% of total) that analytical solution based on 4 part exponential
  // SMM, 17/12/2014
  //============================================================================
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function returns the number of atoms given an effective erosion rate
// (this is the erosion rate in g/cm^2/yr)
// It is used in an optimisation loop (equivalent to fzero in matlab)
// Overloaded, this one returns the spallation and muon contributions
// to the total atoms
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::CRONUS_calculate_N_forward(double effective_erosion_rate, 
                            LSDCRNParameters& LSDCRNP,
                            vector<double>& z_mu, vector<double>& P_mu_z_10Be, 
                            vector<double>& P_mu_z_26Al, double thickSF, 
                            double P_sp_10Be, double P_sp_26Al,
                            double& N_Be10, double& N_Al26, 
                            double& Be10_mu_N, double& Al26_mu_N,
                            double& Be10_sp_N, double& Al26_sp_N)
{                            
  // get the initial fluxes from muons and spallation
  LSDCRNP.integrate_muon_flux_for_erosion(effective_erosion_rate,z_mu, P_mu_z_10Be,
                           P_mu_z_26Al, Be10_mu_N, Al26_mu_N);
                           
  //cout << "Initial muon production guess 10Be: " << Be10_mu_N << " Al26: " << Al26_mu_N << endl;
  LSDCRNP.integrate_nonTD_spallation_flux_for_erosion(effective_erosion_rate,thickSF,
                           P_sp_10Be, P_sp_26Al,Be10_sp_N, Al26_sp_N);

  N_Be10 = Be10_sp_N+Be10_mu_N;
  N_Al26 = Al26_sp_N+Al26_mu_N;
  //cout << "Initial spallation production guess 10Be: " << Be10_sp_N << " Al26: " << Al26_sp_N << endl;                                                                                                                             
  

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function does the error propagation from the CRONUS calculator
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParticle::CRONUS_error_propagation(double pressure, 
                            LSDCRNParameters& LSDCRNP, double thickSF, double rho,
                            double sample_N10, double sample_N26,
                            double sample_del10, double sample_del26,
                            vector<double>& z_mu, vector<double>& P_mu_z_10Be, 
                            vector<double>& P_mu_z_26Al,
                            double P_sp_10Be, double P_sp_26Al, 
                            double eff_e_10, double eff_e_26)
{
  // the vector for holding the uncertainties
  vector<double> uncertainties(12,0.0);

  // get some paramters about production
  double Pmu0_10 = P_mu_z_10Be[0];
  double Pmu0_26 = P_mu_z_26Al[0];
  
  // get the decay coefficients
  bool use_CRONUS = true;
  vector<double> decay_coeff =  LSDCRNP.get_decay_coefficients(use_CRONUS);
  double GammaSp = LSDCRNP.get_spallation_attenuation_length(use_CRONUS);
  
  // get uncertanty parameters for muons
  vector<double> muon_uncertanty_params = 
         LSDCRNP.CRONUS_get_muon_uncertainty_params(pressure);

  // get the relative scalings
  string scaling_name = "St";
  vector<double> rel_delP =
        LSDCRNP.CRONUS_get_uncert_production_ratios(scaling_name);

  //cout << endl << endl << endl << "Entering error checking routine" << endl;
  //cout << "Pmu0_10: " << Pmu0_10 << " Pmu0_26: " << Pmu0_26  << endl;
  //cout << "delPmu_10: " << muon_uncertanty_params[4] << " delPmu_26: " << muon_uncertanty_params[5] << endl;
  //cout << "delPmu_10: " << muon_uncertanty_params[4] << " delPmu_26: " << muon_uncertanty_params[5] << endl;
  //cout << "rel_delP10: " << rel_delP[0] << " rel_delP26: " <<  rel_delP[1] << endl;

  // paramters to calculate individual uncertainties of AMS and production
  double GammaMu_10, GammaMu_26;
  double Be10_mu_N, Al26_mu_N;
  double Psp0_10, Psp0_26;
  double delPsp0_10, delPsp0_26;
  double dEdN_10, dEdN_26;
  double dEdsp0_10, dEdsp0_26;
  double dEdPmu0_10, dEdPmu0_26;
  
  // parameters for error terms
  double delE_ext_10 = 0;
  double delE_int_10 = 0;
  double delE_ext_26 = 0;
  double delE_int_26 = 0;

  // number of atoms determined by the erosion rate function N forward
  double N_Be10, N_Al26, Be10_sp_N, Al26_sp_N;
  
  // variables that are the erosion rate about the centred value
  // used for differencing
  double E_minus, E_plus; 

  // first do 10Be
  if (eff_e_10 > 0)
  {
    // get the contributions from spallation and muons
    CRONUS_calculate_N_forward(eff_e_10, LSDCRNP, z_mu, P_mu_z_10Be, P_mu_z_26Al, 
                               thickSF, P_sp_10Be, P_sp_26Al, N_Be10, N_Al26, 
                               Be10_mu_N, Al26_mu_N, Be10_sp_N, Al26_sp_N);  
     
    uncertainties[8] = Be10_sp_N;
    uncertainties[9] = Be10_mu_N;                         
    
    //cout << "effective e: " << eff_e_10 << " N_Be10: " << N_Be10 
    //     << " N_mu10: " << Be10_mu_N << " N_sp10: " << Be10_sp_N << endl;
                               
    // get the length scale for muons
    GammaMu_10 = eff_e_10/((Pmu0_10/Be10_mu_N)-decay_coeff[0]);
    
    
    // get the production scaling
    Psp0_10 = Be10_sp_N*(decay_coeff[0]+(eff_e_10/GammaSp));
    delPsp0_10 = Psp0_10*rel_delP[0];
    //cout << "GammaMu_10: " << GammaMu_10 << " Psp0_10: " << Psp0_10 
    //     << " delPsp0_10: "<< delPsp0_10 << endl;
    
    // get the values of erosion to do the central difference.
    // This one is for uncertanty in the number of atoms 
    E_plus = CRONUS_simple_N_findroots(eff_e_10, sample_N10+sample_del10,
                              Psp0_10, Pmu0_10, GammaSp, GammaMu_10, decay_coeff[0]);
    E_minus = CRONUS_simple_N_findroots(eff_e_10, sample_N10-sample_del10,
                              Psp0_10, Pmu0_10, GammaSp, GammaMu_10, decay_coeff[0]);
                              
    // NOTE in Greg Balco's code this is multiplied by 10^4 rather than 10^7: the 
    // difference is because Greg's rho is in g/cm^3 and I use kg/m^3
    
    // Greg's code introduces rho here. Presumably this is to convert to m/yr. 
    // WHY???? Everything before has been in g/cm^2 or g/cm^2/yr
    dEdN_10 = (1e7/(2.0*rho*sample_del10))*(E_plus-E_minus);
    
    //cout << "10Be, dN, Eplus: " << E_plus << " E_minus: " << E_minus 
    //     << " dEdn_10: " << dEdN_10 << endl;

    // get the values of erosion to do the central difference. 
    // This one is for the uncertainty in spallation production rate
    E_plus = CRONUS_simple_N_findroots(eff_e_10, sample_N10,
                              Psp0_10+delPsp0_10, Pmu0_10, GammaSp, 
                              GammaMu_10, decay_coeff[0]);
    E_minus = CRONUS_simple_N_findroots(eff_e_10, sample_N10,
                              Psp0_10-delPsp0_10, Pmu0_10, GammaSp, 
                              GammaMu_10, decay_coeff[0]);

    dEdsp0_10 = (1e7/(2.0*rho*delPsp0_10))*(E_plus-E_minus);
    //cout << "10Be, spP0, Eplus: " << E_plus << " E_minus: " << E_minus 
    //     << " dEdps0_10: " << dEdsp0_10 << endl;

    // get the values of erosion to do the central difference. 
    // This one is for the uncertainty in muon production rate
    E_plus = CRONUS_simple_N_findroots(eff_e_10, sample_N10,
                              Psp0_10, Pmu0_10+muon_uncertanty_params[4], GammaSp, 
                              GammaMu_10, decay_coeff[0]);
    E_minus = CRONUS_simple_N_findroots(eff_e_10, sample_N10,
                              Psp0_10, Pmu0_10-muon_uncertanty_params[4], GammaSp, 
                              GammaMu_10, decay_coeff[0]);

    dEdPmu0_10 = (1e7/(2.0*rho*muon_uncertanty_params[4]))*(E_plus-E_minus);
    //cout << "10Be, spP0, Eplus: " << E_plus << " E_minus: " << E_minus 
    //     << " dEdPmu0_10: " << dEdPmu0_10 << endl;
         
         
    delE_ext_10 = sqrt( (dEdsp0_10*delPsp0_10)*(dEdsp0_10*delPsp0_10)
                      + (dEdPmu0_10*muon_uncertanty_params[4])*
                        (dEdPmu0_10*muon_uncertanty_params[4])
                      + (dEdN_10*sample_del10)*(dEdN_10*sample_del10) );
    delE_int_10 = fabs(dEdN_10*sample_del10); 
    
    //cout << "delE_ext_10: " << delE_ext_10 << " delE_int_10: " << delE_int_10 << endl;
    uncertainties[0] = delE_ext_10;
    uncertainties[1] = delE_int_10;
    uncertainties[2] = Psp0_10;
    uncertainties[3] = Pmu0_10;
        
  }
  // then do 26
  if (eff_e_26 > 0)
  {
    // get the contributions from spallation and muons
    CRONUS_calculate_N_forward(eff_e_26, LSDCRNP, z_mu, P_mu_z_10Be, P_mu_z_26Al, 
                               thickSF, P_sp_10Be, P_sp_26Al, N_Be10, N_Al26, 
                               Be10_mu_N, Al26_mu_N, Be10_sp_N, Al26_sp_N);

    // get seperate spallation and muon production for bug checking     
    uncertainties[10] = Al26_sp_N;
    uncertainties[11] = Al26_mu_N; 

    //cout << "effective e: " << eff_e_26 << " N_Al26: " << N_Al26 
    //     << " N_mu26: " << Al26_mu_N << " N_sp26: " << Al26_sp_N << endl;
                               
    // get the length scale for muons
    GammaMu_26 = eff_e_26/((Pmu0_26/Al26_mu_N)-decay_coeff[1]); 
    
    // get the production scaling
    Psp0_26 = Al26_sp_N*(decay_coeff[1]+(eff_e_26/GammaSp));
    delPsp0_26 = Psp0_26*rel_delP[1];
    
    //cout << "GammaMu_26: " << GammaMu_26 << " Psp0_26: " << Psp0_26 
    //     << " delPsp0_26: " << delPsp0_26 << endl;
    
    // get the values of erosion to do the central difference.
    // This one is for uncertanty in the number of atoms  
    E_plus = CRONUS_simple_N_findroots(eff_e_26, sample_N26+sample_del26,
                            Psp0_26, Pmu0_26, GammaSp, GammaMu_26, decay_coeff[1]);
    E_minus = CRONUS_simple_N_findroots(eff_e_26, sample_N26-sample_del26,
                            Psp0_26, Pmu0_26, GammaSp, GammaMu_26, decay_coeff[1]);
                            
    // NOTE in Greg Balco's code this is multiplied by 10^4 rather than 10: the 
    // difference is because Greg's rho is in g/cm^3 and I use kg/m^3
    dEdN_26 = (1.0e7/(2.0*rho*sample_del26))*(E_plus-E_minus);
    
    //cout << "dN, 26AlEplus: " << E_plus << " E_minus: " << E_minus 
    //     << " dEdn_10: " << dEdN_26 << endl;

    // get the values of erosion to do the central difference. 
    // This one is for the uncertainty in the production rate
    E_plus = CRONUS_simple_N_findroots(eff_e_26, sample_N26,
                              Psp0_26+delPsp0_26, Pmu0_26, GammaSp, 
                              GammaMu_26, decay_coeff[1]);
    E_minus = CRONUS_simple_N_findroots(eff_e_26, sample_N26,
                              Psp0_26-delPsp0_26, Pmu0_26, GammaSp, 
                              GammaMu_26, decay_coeff[1]);
    
    dEdsp0_26 = (1e7/(2.0*rho*delPsp0_26))*(E_plus-E_minus);
    //cout << "26Al, spP0, Eplus: " << E_plus << " E_minus: " << E_minus 
    //     << " dEdsp_26: " << dEdsp0_26 << endl;

    // get the values of erosion to do the central difference. 
    // This one is for the uncertainty in muon production rate
    E_plus = CRONUS_simple_N_findroots(eff_e_26, sample_N26,
                              Psp0_26, Pmu0_26+muon_uncertanty_params[5], GammaSp, 
                              GammaMu_26, decay_coeff[1]);
    E_minus = CRONUS_simple_N_findroots(eff_e_26, sample_N26,
                              Psp0_26, Pmu0_26-muon_uncertanty_params[5], GammaSp, 
                              GammaMu_26, decay_coeff[1]);

    dEdPmu0_26 = (1e7/(2.0*rho*muon_uncertanty_params[5]))*(E_plus-E_minus);
    //cout << "26Be, spP0, Eplus: " << E_plus << " E_minus: " << E_minus 
    //     << " dEdPmu0_26: " << dEdPmu0_26 << endl;
         
    delE_ext_26 = sqrt( (dEdsp0_26*delPsp0_26)*(dEdsp0_26*delPsp0_26)
                      + (dEdPmu0_26*muon_uncertanty_params[5])*
                        (dEdPmu0_26*muon_uncertanty_params[5])
                      + (dEdN_26*sample_del26)*(dEdN_26*sample_del26) );
    delE_int_26 = fabs(dEdN_26*sample_del26);           
    //cout << "delE_ext_26: " << delE_ext_26 << " delE_int_26: " << delE_int_26 << endl;   

    uncertainties[4] = delE_ext_26;
    uncertainties[5] = delE_int_26; 
    uncertainties[6] = Psp0_26;
    uncertainties[7] = Pmu0_26;
  }  
  
  return uncertainties;



}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis finds the root of the simple erosion rate finder, used to 
// Replicate CRONUS uncertanty analysis
// The function uses Newton-Raphs, so may give slightly different answers than
// the CRONUS version, which uses Matlab's fzero function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::CRONUS_simple_N_findroots(double inital_erate, double target,
                                           double Psp, double Pmu, double GammaSp,  
                                           double GammaMu,double decay)
{
  // paramters controlling the erosion rates
  double e_new;
  double this_erate;
  e_new = inital_erate;
  double e_displace = 1e-5;
  double displace_erate;
  double e_change;
  
  // parameters for Newton-Raphson
  double f_x;
  double f_x_displace;
  double e_derivative;
  double tolerance = 1e-10;

  // now iterate using newton-raphson
  do
  {
    
    this_erate = CRONUS_simple_N_forward(e_new, Psp, Pmu, GammaSp, GammaMu, decay);

     f_x =  this_erate-target;                  
                        
    // now get the derivative
    displace_erate = CRONUS_simple_N_forward(e_new+e_displace, Psp, Pmu, GammaSp, 
                                             GammaMu, decay);
    f_x_displace =  displace_erate-target;                    
                        
    e_derivative = (f_x_displace-f_x)/e_displace;
    
    //cout << "N10: " << N10_this_step << " + displace: " << N10_displace  << endl;
    //cout << "fx: " <<  f_x << " f_x_displace: " << f_x_displace << " fprimex: " << N_derivative << endl;
    
    if(e_derivative != 0)
    {
      e_new = e_new-f_x/e_derivative;
    
      // check to see if the difference in erosion rates meet a tolerance
      e_change = f_x/e_derivative;
      //cout << "Change is: " << e_change << " and erosion rate is: " << e_new << endl;
    }
    else
    {
      e_change = 0;
    }
  
  } while(fabs(e_change) > tolerance);
  return e_new; 

}                             
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Simple N forward: a function for calculating erosion that is used in the 
// error propagation routines
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::CRONUS_simple_N_forward(double eff_eros, double Psp, double Pmu,
                                      double GammaSp, double GammaMu, double decay)
{
  double N = (Psp/(decay + eff_eros/GammaSp)) 
              + (Pmu/(decay + eff_eros/GammaMu));
  return N;            

}                                      
#endif


