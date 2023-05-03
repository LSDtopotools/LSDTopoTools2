///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDCosmoRaster.cpp
/// header for the CosmoRaster object
/// The cosmo raster contains routines for calculating cosmogenic concentrations
/// from a variety of erosion scenarios
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1    17/11/2018
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDCRNParameters.hpp"
using namespace std;
using namespace TNT;

// Sorting compiling problems with MSVC
#ifdef _WIN32
#ifndef M_PI
extern double M_PI;
#endif
#endif

#ifndef LSDCosmoRaster_H
#define LSDCosmoRaster_H

///@brief Create model objects to use LSDRaster methods on synthetic landscapes.
class LSDCosmoRaster: public LSDRaster
{
  public:

    /// @brief Constructor. Create an LSDCosmoRaster from an LSDRaster.
    /// @return LSDCosmoRaster
    LSDCosmoRaster()
    {
      create();
    }


    /// @brief Constructor. Create an LSDCosmoRaster from an LSDRaster.
    /// @return LSDCosmoRaster
    /// @param An_LSDRaster LSDRaster object.
    LSDCosmoRaster(LSDRaster& An_LSDRaster)
    {
      create(An_LSDRaster);
    }



    /// @brief Constructor. Create an LSDCosmoRaster from a file.
    /// Uses a filename and file extension
    /// @return LSDCosmoRaster
    /// @param filename A String, the file to be loaded.
    /// @param extension A String, the file extension to be loaded.
    LSDCosmoRaster(string filename, string extension)
    {
      create(filename, extension);
    }

    /// @brief Constructor. Create a blank raster nodel
    /// @return LSDRasterModel
    /// @param NCols Height of raster
    /// @param NRows Width of raster
    LSDCosmoRaster(int NRows, int NCols)
    {
      create(NRows,NCols);
    }

    /// @brief This just returns the raster model object data as a raster
    /// @return A raster with the data from the LSDRasterModel
    /// @author SMM
    /// @date 01/09/2017
    LSDRaster return_as_raster();

    /// @brief This resizes the LSDRasterModel, resetting some flags in the process,
    /// as well as setting many of the Array2D data members to be empty arrays
    /// The raster data in the end is a random surface (determined by the noise
    /// data member)
    /// This overloaded version also changes the data resolution
    /// @param new_rows the new number of rows
    /// @param new_cols the new number of columns
    /// @param new_resolution the new data resolution
    /// @param new_value the new value of the raster
    /// @author SMM
    /// @date 01/09/2017
    void resize_and_reset( int new_rows, int new_cols, float new_resolution, float new_value );

    /// @brief Gets the minimum and maximum values from a raster
    /// @return vector with the first element is the minimum and second is the maximum
    /// @author SMM
    /// @date 03/09/2017
    vector<float> minimum_and_maximum_value();

    /// @brief This returns a clipped raster that has the same dimensions as the
    ///  smaller raster
    /// @param smaller_raster the raster to which the bigger raster should be
    ///  clipped
    /// @author SMM
    /// @date 20/03/2015
    LSDRaster clip_to_smaller_raster(LSDRaster& smaller_raster);

    /// @brief This returns a clipped raster that has the same dimensions as the
    ///  smaller raster
    /// @param smaller_raster the raster to which the bigger raster should be
    ///  clipped
    /// @author SMM
    /// @date 20/03/2015
    LSDRaster clip_to_smaller_raster(LSDIndexRaster& smaller_raster);

    ///@brief This function returns the raster data as text file
    ///@return text file with raster data
    ///@author FJC
    ///@date 30/09/16
    void write_RasterData_to_text_file(string filename);

    ///@brief This function returns the raster data as a vector
    ///@return vector<float> with raster data
    ///@author FJC
    ///@date 06/11/15
    vector<float> get_RasterData_vector();

    /// @brief This function returns a vector with the X adn Y minimum and max
    ///   values
    /// @return XYMinMax a vector with four elements
    ///   XYMinMax[0] = XMinimum
    ///   XYMinMax[1] = YMinimum
    ///   XYMinMax[2] = XMaximum
    ///   XYMinMax[3] = XMaximum
    /// @author SMM
    /// @date 3/7/2015
    vector<float> get_XY_MinMax();
  
    /// @brief This resets scaling of the CRN parameters object
    /// @param LSDCRNP an LSDCRNParameters object. Is passed by refernce and reset in the function
    /// @param muon scaling. A string. Options are: Schaller, Braucher, Granger, newCRONUS.
    ///    the default is Braucher but it will spit out a warning if you don't assign this. 
    /// @author SMM
    /// @date 18/11/2018
    void reset_scaling(LSDCRNParameters& LSDCRNP, string Muon_scaling);
  
    /// @brief This creates a production raster from elevation data that is 
    ///  later used in scaling calculations
    /// @return a raster with production scaling
    /// @param Elevation_Data an elevation raster
    /// @param path_to_atmospheric_data A strong containing the path tothe atmospheric NCEP data
    /// @author SMM
    /// @date 18/11/2018
    LSDRaster calculate_production_raster(LSDRaster& Elevation_Data, string path_to_atmospheric_data);

    /// @brief This creates an atmospheric pressure raster for debugging 
    /// @return a raster with atmospheric pressure
    /// @param Elevation_Data an elevation raster
    /// @param path_to_atmospheric_data A strong containing the path tothe atmospheric NCEP data
    /// @author SMM
    /// @date 18/11/2018
    LSDRaster calculate_atm_pressure_raster(LSDRaster& Elevation_Data, string path_to_atmospheric_data);


    /// @brief this predicts the concentration of a nuclide for each pixel within a basin. 
    /// It does a full analytical solution to account for a step change in the erosion rate
    /// there are functions here that calculate timing of step change so use those to prepare 
    /// input rasters
    /// @detail This is NOT an accumulator function! It is the concentration at that specific pixel. 
    ///  For pixels with a landslide (e.g, places with TopoShield not = 0)
    ///  This is the depth averaged concentration of particles emerging from that pixel during the landslide. 
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger, and newCRONUS
    /// @param erosion_rate_m_yr A raster of the erosion rate (in m/yr)
    /// @param time_since_change The time since the change in erosion rate
    /// @param ProductionScale A raster of the production scaling (uses Lal/Stone scheme)
    /// @param TopoShield A raster of the topographic shielding. Varies from 0 to 1.
    /// @param old_erosion_rate_m_yr Erosion rate before step change in m/yr
    /// @param new_erosion_rate_m_yr Erosion rate after step change in m/yr  
    /// @param rock_density_kg_m3 rock density to convert to g/cm^2/yr for the effective erosion rates
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide at each pixel for the given erosion rate
    /// @author SMM
    /// @date 03/06/2022
    LSDRaster calculate_CRN_concentration_raster_step_change(string Nuclide,
                                           string Muon_scaling, 
                                           LSDRaster& erosion_rate_m_yr,
                                           LSDRaster& time_since_change,
                                           LSDRaster& ProductionScale, 
                                           LSDRaster& TopoShield,  
                                           float old_erosion_rate_m_yr, float new_erosion_rate_m_yr,
                                           float rock_density_kg_m3,
                                           bool is_production_uncertainty_plus_on,
                                           bool is_production_uncertainty_minus_on);


    /// @brief this predicts the concentration of a nuclide for each pixel within a basin. 
    /// It does a full analytical solution to account for
    ///  snow and self sheilding
    /// @detail This is NOT an accumulator function! It is the concentration at that specific pixel. 
    ///  For pixels with a landslide (e.g, places with TopoShield not = 0)
    ///  This is the depth averaged concentration of particles emerging from that pixel during the landslide. 
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger, and newCRONUS
    /// @param eff_erosion_rate A raster of the effective erosion rate (in g/cm^2/yr)
    /// @param ProductionScale A raster of the production scaling (uses Lal/Stone scheme)
    /// @param TopoShield A raster of the topographic shielding. Varies from 0 to 1.
    /// @param SelfShield A raster of self shielding in g/cm^2 (the effective thickenss of the layer removed) 
    /// @param SnowShield A raster of snow shielding in g/cm^2 (the effective depth of the snow) 
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide at each pixel for the given erosion rate
    /// @author SMM
    /// @date 18/11/2018
    LSDRaster calculate_CRN_concentration_raster(string Nuclide, string Muon_scaling, LSDRaster& eff_erosion_rate,
                               LSDRaster& ProductionScale, LSDRaster& TopoShield, 
                               LSDRaster& SelfShield, LSDRaster& SnowShield, 
                               bool is_production_uncertainty_plus_on,
                               bool is_production_uncertainty_minus_on);

    /// @brief this predicts the concentration of a nuclide for each pixel within a basin. 
    /// It does a full analytical solution to account for
    ///  snow and self sheilding. This is overloaded and only calculates the concentrations
    ///  based on an outlet node
    /// @detail This is NOT an accumulator function! It is the concentration at that specific pixel. 
    ///  For pixels with a landslide (e.g, places with TopoShield not = 0)
    ///  This is the depth averaged concentration of particles emerging from that pixel during the landslide. 
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger, and newCRONUS
    /// @param eff_erosion_rate A raster of the effective erosion rate (in g/cm^2/yr)
    /// @param ProductionScale A raster of the production scaling (uses Lal/Stone scheme)
    /// @param TopoShield A raster of the topographic shielding. Varies from 0 to 1.
    /// @param SelfShield A raster of self shielding in g/cm^2 (the effective thickenss of the layer removed) 
    /// @param SnowShield A raster of snow shielding in g/cm^2 (the effective depth of the snow) 
    /// @param FlowInfo A flow info object
    /// @param upslope nodes The nodes upslope of an arbitrary outlet node. 
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide at each pixel for the given erosion rate
    /// @author SMM
    /// @date 31/07/2019
    LSDRaster calculate_CRN_concentration_raster(string Nuclide, string Muon_scaling, LSDRaster& eff_erosion_rate,
                               LSDRaster& ProductionScale, LSDRaster& TopoShield, 
                               LSDRaster& SelfShield, LSDRaster& SnowShield, 
                               LSDFlowInfo& FlowInfo, int outlet_node,
                               bool is_production_uncertainty_plus_on,
                               bool is_production_uncertainty_minus_on);


    /// @brief this predicts the concentration of a nuclide for each pixel within a basin. 
    /// It does a full analytical solution to account for
    ///  snow and self sheilding. This is overloaded and only calculates the concentrations
    ///  based on a node list derived from a flowinfo object. 
    /// @detail This is NOT an accumulator function! It is the concentration at that specific pixel. 
    ///  For pixels with a landslide (e.g, places with TopoShield not = 0)
    ///  This is the depth averaged concentration of particles emerging from that pixel during the landslide. 
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger, and newCRONUS
    /// @param eff_erosion_rate A raster of the effective erosion rate (in g/cm^2/yr)
    /// @param ProductionScale A raster of the production scaling (uses Lal/Stone scheme)
    /// @param TopoShield A raster of the topographic shielding. Varies from 0 to 1.
    /// @param SelfShield A raster of self shielding in g/cm^2 (the effective thickenss of the layer removed) 
    /// @param SnowShield A raster of snow shielding in g/cm^2 (the effective depth of the snow) 
    /// @param FlowInfo A flow info object
    /// @param upslope nodes The nodes upslope of an arbitrary outlet node. 
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide at each pixel for the given erosion rate
    /// @author SMM
    /// @date 29/07/2019
    LSDRaster calculate_CRN_concentration_raster(string Nuclide, string Muon_scaling, LSDRaster& eff_erosion_rate,
                               LSDRaster& ProductionScale, LSDRaster& TopoShield, 
                               LSDRaster& SelfShield, LSDRaster& SnowShield, 
                               LSDFlowInfo& FlowInfo, vector<int> upslope_nodes,
                               bool is_production_uncertainty_plus_on,
                               bool is_production_uncertainty_minus_on);

    /// @brief This returns a raster that has the erosion rate including a "dump" of sediment from landslides
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param SelfShield A raster of self shielding in g/cm^2 (the effective thickenss of the layer removed) 
    /// @return The rate (in g/cm^2/yr) where the "dumped" sediment appears in the rate (e.g., assumed to be delived to the system over a year)
    ///  or perhaps a better way to think about it is that the sediment in the system is at the background erosion rate and one year's worth of this
    ///  sediment is diluted by the landslide material. 
    /// @author SMM
    /// @date 30/07/2019    
    LSDRaster calculate_landslide_dump_erate(LSDRaster& eff_erosion_rate,
                                           LSDRaster& SelfShield);

    /// @brief This accumulates CRN and gets detrital CRN concentrations. It is a heavily overloaded function
    ///  This version assumes all quartz content is the same and does every pixel in the DEM
    /// @param CRN_conc The concentration of CRN derived from the 
    ///   calculate_CRN_concentration_raster function
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param FlowInfo A flow info object
    /// @return the concentration of the nuclide at each pixel for the given erosion rate that is accumulated in detrital sediment
    /// @author SMM
    /// @date 30/07/2019
    LSDRaster calculate_accumulated_CRN_concentration(LSDRaster& CRN_conc, 
                                          LSDRaster& eff_erosion_rate,
                                          LSDFlowInfo& FlowInfo);

    /// @brief This accumulates CRN and gets detrital CRN concentrations. It is a heavily overloaded function
    ///  This version takes a quartz concentration raster
    /// @param CRN_conc The concentration of CRN derived from the 
    ///   calculate_CRN_concentration_raster function
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param quartz_concentration The quartz "concentration": in fact this is the relative mass fraction of quartz in eroded material 
    ///  in each pixel. 
    /// @param FlowInfo A flow info object
    /// @return the concentration of the nuclide at each pixel for the given erosion rate that is accumulated in detrital sediment
    /// @author SMM
    /// @date 30/07/2019
    LSDRaster calculate_accumulated_CRN_concentration(LSDRaster& CRN_conc, 
                                          LSDRaster& eff_erosion_rate,
                                          LSDRaster& quartz_concentration,
                                          LSDFlowInfo& FlowInfo);                                      

    /// @brief This accumulates CRN and gets detrital CRN concentrations. It is a heavily overloaded function
    ///  This is for a single pixel
    /// @param CRN_conc The concentration of CRN derived from the 
    ///   calculate_CRN_concentration_raster function
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param quartz_concentration The quartz "concentration": in fact this is the relative mass fraction of quartz in eroded material 
    ///  in each pixel. 
    /// @param FlowInfo A flow info object
    /// @param upslope_nodes_and_outlet node indices of uplsope nodes and the outlet. 
    ///  get from LSDFlowInfo::get_uplsope_nodes_include_outlet
    /// @return the concentration of the accumulated nuclide at the outlet node
    /// @author SMM
    /// @date 30/07/2019
    float calculate_accumulated_CRN_concentration(LSDRaster& CRN_conc, 
                                          LSDRaster& eff_erosion_rate,
                                          LSDFlowInfo& FlowInfo, 
                                          vector<int> upslope_nodes_and_outlet);    

    /// @brief This accumulates CRN and gets detrital CRN concentrations. It is a heavily overloaded function
    ///  This version takes a quartz concentration raster and is for a single pixel
    /// @param CRN_conc The concentration of CRN derived from the 
    ///   calculate_CRN_concentration_raster function
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param quartz_concentration The quartz "concentration": in fact this is the relative mass fraction of quartz in eroded material 
    ///  in each pixel. 
    /// @param FlowInfo A flow info object
    /// @param upslope_nodes_and_outlet node indices of uplsope nodes and the outlet. 
    ///  get from LSDFlowInfo::get_uplsope_nodes_include_outlet
    /// @return the concentration of the accumulated nuclide at the outlet node
    /// @author SMM
    /// @date 30/07/2019
    float calculate_accumulated_CRN_concentration(LSDRaster& CRN_conc, 
                                          LSDRaster& eff_erosion_rate,
                                          LSDRaster& quartz_concentration,
                                          LSDFlowInfo& FlowInfo, 
                                          vector<int> upslope_nodes_and_outlet);    

    /// @brief This accumulates CRN and gets detrital CRN concentrations. It is a heavily overloaded function
    ///  This version is for a single pixel
    /// @param CRN_conc The concentration of CRN derived from the 
    ///   calculate_CRN_concentration_raster function
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param quartz_concentration The quartz "concentration": in fact this is the relative mass fraction of quartz in eroded material 
    ///  in each pixel. 
    /// @param FlowInfo A flow info object
    /// @param outlet_node node index of outlet node. 
    /// @return the concentration of the accumulated nuclide at the outlet node
    /// @author SMM
    /// @date 30/07/2019
    float calculate_accumulated_CRN_concentration(LSDRaster& CRN_conc, 
                                          LSDRaster& eff_erosion_rate,
                                          LSDFlowInfo& FlowInfo, 
                                          int outlet_node);    

    /// @brief This accumulates CRN and gets detrital CRN concentrations. It is a heavily overloaded function
    ///  This version takes a quartz concentration raster and is for a single pixel
    /// @param CRN_conc The concentration of CRN derived from the 
    ///   calculate_CRN_concentration_raster function
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param quartz_concentration The quartz "concentration": in fact this is the relative mass fraction of quartz in eroded material 
    ///  in each pixel. 
    /// @param FlowInfo A flow info object
    /// @param outlet_node node index of outlet node. 
    /// @return the concentration of the accumulated nuclide at the outlet node
    /// @author SMM
    /// @date 30/07/2019
    float calculate_accumulated_CRN_concentration(LSDRaster& CRN_conc, 
                                          LSDRaster& eff_erosion_rate,
                                          LSDRaster& quartz_concentration,
                                          LSDFlowInfo& FlowInfo, 
                                          int outlet_node);                                         
 
    /// @brief this predicts the mean concentration of a nuclide for each pixel within a basin. 
    /// It does a full analytical solution to account for
    ///  snow and self sheilding
    /// @detail This is NOT an accumulator function! It is the concentration at that specific pixel. 
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger, and newCRONUS
    /// @param eff_erosion_rate A raster of the effective erosion rate (in g/cm^2/yr)
    /// @param ProductionScale A raster of the production scaling (uses Lal/Stone scheme)
    /// @param TopoShield A raster of the topographic shielding. Varies from 0 to 1.
    /// @param SelfShield A raster of self shielding in g/cm^2 (the effective thickenss of the layer removed) 
    /// @param SnowShield A raster of snow shielding in g/cm^2 (the effective depth of the snow) 
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 18/11/2018
    void predict_mean_CRN_conc(string Nuclide, string Muon_scaling, LSDRaster& eff_erosion_rate,
                               LSDRaster& ProductionScale, LSDRaster& TopoShield, 
                               LSDRaster& SelfShield, LSDRaster& SnowShield, 
                               bool is_production_uncertainty_plus_on,
                               bool is_production_uncertainty_minus_on);


                                       

    /// @brief This uses newton-raphson iteration to get the erosion rate in unknown pixels
    ///  that reproduced a known nuclide concentration
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger, and newCRONUS
    /// @param known_eff_erosion_rate A raster of the effective erosion rate (in g/cm^2/yr) where it is known
    ///  nodata elsewhere
    /// @param ProductionScale A raster of the production scaling (uses Lal/Stone scheme)
    /// @param TopoShield A raster of the topographic shielding. Varies from 0 to 1.
    /// @param SelfShield A raster of self shielding in g/cm^2 (the effective thickenss of the layer removed) 
    /// @param SnowShield A raster of snow shielding in g/cm^2 (the effective depth of the snow) 
    /// @paramquartz_concentration The relative concentration of quartz
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @param predicted_nuclide_conc this is return by reference so it takes
    ///   it is replaced by the function. The predicted cosmo concentration
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 18/11/2018
    float calculate_eff_erate_from_conc(float Nuclide_conc,
                                    string Nuclide,string Muon_scaling, 
                                    LSDRaster& known_eff_erosion_rate,
                                    LSDRaster& ProductionScale,
                                    LSDRaster& TopoShield, 
                                    LSDRaster& SelfShield,
                                    LSDRaster& SnowShield,
                                    LSDRaster& quartz_concentration,
                                    LSDFlowInfo& FlowInfo,
                                    int outlet_node,
                                    float& predicted_nuclide_conc);


    /// @brief This uses newton-raphson iteration to get the erosion rate in unknown pixels
    ///  that reproduced a known nuclide concentration
    /// @param Nuclide_conc the concentration (in atoms/g) of the nuclide
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger, and newCRONUS
    /// @param Production_raster A raster of the production scaling (uses Lal/Stone scheme)
    /// @param TopoShield A raster of the topographic shielding. Varies from 0 to 1.
    /// @param SelfShield A raster of self shielding in g/cm^2 (the effective thickenss of the layer removed) 
    /// @param SnowShield A raster of snow shielding in g/cm^2 (the effective depth of the snow) 
    /// @param FlowInfo A flow info object
    /// @param outlet_node node index of outlet node. 
    /// @return an effective erosion rate from the outlet pixel at production rates of that pixel
    /// @author SMM
    /// @date 09/08/2019
    float get_erate_guess(float Nuclide_conc, string Nuclide,string Muon_scaling, 
                                      LSDRaster& Production_raster,
                                      LSDRaster& TopoShield, 
                                      LSDRaster& SelfShield,
                                      LSDRaster& SnowShield,
                                      LSDFlowInfo& FlowInfo,
                                      int outlet_node);


    /// @brief This creates an erate raster that has a fixed erate everywhere except where
    ///  there is a known erate (calculated, for example, using a nested cosmo sample)
    /// @param eff_erate The e4ffective erate that is used in places where there are no "known" erates (in g/cm^2/yr)
    /// @param known_eff_erosion_rate A raster of the effective erosion rate (in g/cm^2/yr) where it is known
    ///  nodata elsewhere
    /// @param outlet_node node index of outlet node. 
    /// @return The raster with effective erates
    /// @author SMM
    /// @date 09/08/2019
    LSDRaster make_eff_erate_raster_from_known_erate_and_float_erate(float eff_erate, 
                                    LSDRaster& known_erate_raster,
                                    int outlet_node,
                                    LSDFlowInfo& FlowInfo);                                   


    /// @brief This returns the rasters for a step change 
    /// @author SMM
    /// @date 02/06/2022
    void calculate_step_change_rasters(LSDFlowInfo& FlowInfo, vector<int>& basin_nodes, LSDRaster& chi_coordinate, LSDRaster& erate, LSDRaster& set_topography, 
                                        LSDRaster& time_since_step, float time_since_start_of_step, float old_uplift, float new_uplift,
                                        float chi_knickpoint, float z_knickpoint, float n, float ksn_old, float ksn_new, 
                                        float z_outlet, float K);  

    /// @brief This returns the rasters for a step change 
    /// @author SMM
    /// @date 03/06/2022
    void calculate_step_change_hillslopes(LSDRaster& fluvial_topography, LSDRaster& raster_tags, LSDRaster& new_topography, 
                                                      LSDRaster& erate, LSDRaster& time_since_step, vector<string> boundary_conditions, 
                                                      int threshold_pixels, float erate_old, float Sc_value_old, float Sc_value_new);                                      


  protected:



  private:
    /// @brief Make a 100x100 raster
    void create();

    /// @brief Load a raster
    void create(string filename, string extension);

    /// @brief Make a CosmoRaster from another raster
    void create(LSDRaster& An_LSDRaster);

    /// @brief create a raster from NRows and NCols information
    void create(int NRows, int NCols);

};

#endif
