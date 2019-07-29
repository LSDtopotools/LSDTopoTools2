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
  
    /// @ brief This resets scaling of the CRN parameters object
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
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 18/11/2018
    void calculate_CRN_concentration_raster(string Nuclide, string Muon_scaling, LSDRaster& eff_erosion_rate,
                               LSDRaster& ProductionScale, LSDRaster& TopoShield, 
                               LSDRaster& SelfShield, LSDRaster& SnowShield, 
                               bool is_production_uncertainty_plus_on,
                               bool is_production_uncertainty_minus_on);


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
