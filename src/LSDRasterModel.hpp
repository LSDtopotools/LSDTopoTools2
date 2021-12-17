///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDRasterModel.cpp
/// cpp file for the LSDRasterModel object
/// LSD stands for Land Surface Dynamics
/// This object provides an environment for landscape evolution modelling, which can then
/// be integrated with the topographic analysis tools to efficiently analyse model runs.
///
/// The landscape evolution model uses implicit methods to provide stability with
/// relatively long timesteps.  Fluvial erosion is solved following Braun and Willet (2013)
/// using the fastscape algorithm, whilst hillslope sediment transport is modelled as a
/// non-linear diffusive sediment flux, following the implicit scheme developed for
/// MuDDPile.
///
/// The aim is to have two complimentary models:
/// i) a simple coupled hillslope-channel model in which large scale landscape dynamics
/// can be modelled
/// ii) a more complex treatment of hillslopes explicitly incorporating the role of
/// vegetation in driving sediment production and transport, and that copes with the
/// with the transition from soil mantled-bedrock hillslopes at high erosion rates.
///
/// In order to run the model, one needs a parameter file that should be read by the
/// driver function.
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
/// @author David T. Milodowski, University of Edinburgh
/// @author Martin D. Hurst, British Geological Survey
/// @author Fiona Clubb, University of Edinburgh
/// @author Stuart Grieve, University of Edinburgh
/// @author James Jenkinson, University of Edinburgh
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1    24/07/2013
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TNT/tnt.h"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "LSDRaster.hpp"
#include "LSDRasterSpectral.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDSpatialCSVReader.hpp"
#include "LSDParticleColumn.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDLithoCube.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRasterModel_H
#define LSDRasterModel_H

///@brief Create model objects to use LSDRaster methods on synthetic landscapes.
class LSDRasterModel: public LSDRasterSpectral
{
  public:

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  // CONSTRUCTORS AND CREATE FUNCTIONS
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief Constructor. Create a deafult LSDRasterModel (100x100)
  /// @return An instance of LSDRasterModel
  LSDRasterModel()
  {
    create();
  }

  /// @brief Constructor. Create a LSDRasterModel from a parameter file
  /// @return An instance of LSDRasterModel
  /// @param master_param A filenam for the master parameter file
  LSDRasterModel( string master_param )
  {
    create(master_param);
  }

  /// @brief Constructor. Create an LSDRasterModel from a file.
  /// Uses a filename and file extension
  /// @return LSDRasterModel
  /// @param filename A String, the file to be loaded.
  /// @param extension A String, the file extension to be loaded.
  LSDRasterModel(string filename, string extension)
  {
    create(filename, extension);
    default_parameters();
  }

  /// @brief Constructor. Create an LSDRasterModel from memory.
  /// @return LSDRasterModel
  /// @param nrows An integer of the number of rows.
  /// @param ncols An integer of the number of columns.
  /// @param xmin A float of the minimum X coordinate.
  /// @param ymin A float of the minimum Y coordinate.
  /// @param cellsize A float of the cellsize.
  /// @param ndv An integer of the no data value.
  /// @param data An Array2D of floats in the shape nrows*ncols,
  ///containing the data to be written.
  LSDRasterModel(int nrows, int ncols, float xmin, float ymin,
            float cellsize, float ndv, Array2D<float> data, map<string,string> GRS)
  {
    default_parameters();
    create(nrows, ncols, xmin, ymin, cellsize, ndv, data, GRS);
  }

  /// @brief Constructor. Create an LSDRasterModel from an LSDRaster.
  /// @return LSDRasterModel
  /// @param An_LSDRaster LSDRaster object.
  LSDRasterModel(LSDRaster& An_LSDRaster)
  {
    create(An_LSDRaster);
    default_parameters();
  }

  /// @brief Constructor. Create a blank raster nodel
  /// @return LSDRasterModel
  /// @param NCols Height of raster
  /// @param NRows Width of raster
  LSDRasterModel(int NRows, int NCols);

  /// @brief Class destructor
  ~LSDRasterModel( void );

  /// @brief operator assignment
  LSDRasterModel& operator=(const LSDRasterModel& LSDR);

  /// @brief This just returns the raster model object data as a raster
  /// @return A raster with the data from the LSDRasterModel
  /// @author SMM
  /// @date 01/09/2017
  LSDRaster return_as_raster();

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  // INITIALISATION ROUTINES
  // @~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief this initialises the model by directly #ting the data members
  void initialize_model(
    string& parameter_file, string& run_name, float& dt, float& EndTime, float& PrintInterval,
    float& k_w, float& b, float& m, float& n, float& K, float& ErosionThreshold,
    float& K_nl, float& S_c, float& UpliftRate, float& PrecipitationRate,
    float& NorthBoundaryElevation, float& SouthBoundaryElevation,
    Array2D<float>& PrecipitationFlux, Array2D<float>& SlopesBetweenRows,
    Array2D<float>& SlopesBetweenColumns, Array2D<float>& ErosionRate);

  /// @brief This module initialises the model runs, calling the required function from
  /// the initial topography and loads the parameters from the parameter file.
  /// @param parameter_file the filename of the paramter file (with extension)
  /// @author JAJ
  /// @date 01/01/2014
  void initialize_model( string parameter_file );

  /// @brief This module initialises the model using the maps that have been read
  /// in from the parameter file via the param file parser function.
  /// @param nowt
  /// @author DAV
  /// @date 2015-01-17
  void initialise_model();

  /// @brief this appends a string to the run name
  /// @detail can be used to append parameters to run names
  /// @author SMM
  /// @date 09/04/2015
  void append_run_name(string append_name);

  /// @brief Adds random noise to each pixel in range [min, max]
  /// @param minium random addition
  /// @param maximum random addition
  /// @author JAJ
  /// @date 01/01/2014
  void random_surface_noise( float min, float max );

  /// @brief Adds random noise to each pixel using the noise data member
  /// checks no data and you feed it a see
  /// @param seed the seed to the random number (if you give the same seed you should get the same results)
  /// @author SMM
  /// @date 17/06/2014
  void random_surface_noise(long seed);

  /// @brief Adds random noise to each pixel using the noise data member
  /// @author SMM
  /// @date 17/06/2014
  void random_surface_noise();

  /// @brief This resets the RasterData array to have a parabolic shape,
  /// with 0 elevation at the N and S boundaries. It also adds some random
  /// noise to the topography. The amplitude of this noise is set by the
  /// data member 'noise'. The default noise is 1mm.
  /// @param peak_elev The peak elevation in metres. Is in the middle of the
  /// model domain
  /// @param edge_offset an offest from the edge elevation. You can have a little
  /// cliff at the edge
  /// @author SMM
  /// @date 1/7/2014
  void initialise_parabolic_surface(float peak_elev, float edge_offset);

  /// @brief Adds a parabolic surface to the DEM. Used to try and avoid
  ///  ;arge areas of fill from the fractal initiation steps
  /// @param peak_elev The peak elevation in metres. Is in the middle of the
  /// model domain
  /// @author SMM
  /// @date 11/8/2017
  void superimpose_parabolic_surface(float peak_elev);


  /// @brief This initialises the raster model to a square model domain
  /// with a fractal surface using the algorithm from Saupe (1987d)
  /// @param fractal_D Used to determine the fractal dimension, D by:
  /// D = 3 - fractal_D. So the fractal dimension should be between
  /// 2 and 3.
  /// @author DAV
  /// @date 20/10/2014
  void intialise_fourier_fractal_surface(float fractal_D);

  /// @brief This initialises the raster model to a square model domain
  /// with a fractal surface using the algorithm from the LSDRasterModel
  /// @param beta Used to determine the fractal dimension, beta by:
  /// beta = 3 - fractal_beta. So the fractal dimension should be between
  /// 2 and 3.
  /// @param desired_relief The relief desired from the final surface
  /// @author SMM
  /// @date 10/08/2017
  void intialise_fourier_fractal_surface_v2(float beta, float desired_relief);

  /// @brief This initialises the raster model to a fractal surface using the
  ///  diamond square algorithm
  /// @param feature_order is an interger n where the feature size consists of 2^n nodes.
  /// If the feature order is set bigger than the dimensions of the parent raster then
  /// this will default to the order of the parent raster.
  /// @param desired_relief The relief desired from the final surface
  /// @author SMM
  /// @date 10/08/2017
  void intialise_diamond_square_fractal_surface(int feature_order, float desired_relief);

  /// @brief This takes a raster and tapers the edges to zero elevation. Its
  ///  purpose is to remove edge artefacts.
  /// @param rows_to_taper The number of rows at the N and S boudaries to taper
  /// @author SMM
  /// @date 10/08/2017
  void initialise_taper_edges_and_raise_raster(int rows_to_taper);

  /// @brief this raises the raster so the lowest point is zero and also fills the raster
  /// @author SMM
  /// @date 25/08/2017
  void raise_and_fill_raster();
  void raise_and_fill_raster(float min_slope_for_fill);
  void fill_raster(float min_slope_for_fill);

  /// @brief CHanges the elevation of all nodata nodes by adjustment
  /// @param adjustment the elevation to change
  /// @author SMM
  /// @date 08/02/2021
  void add_fixed_elevation(float adjustment);
  
  /// @brief This takes the surface topography and subtracts from it the cumulative uplift, 
  ///  in order to calculate the total exhumation and surface for the lithocube
  /// @param Cumulative_uplift the raster containing the cumulative uplift of this step in m
  /// @return a raster with the exhumation surface for the lithocube
  /// @date 06/05/2021
  /// @author SMM
  LSDRaster calculate_exhumation_surface_from_cumulative_uplft(LSDRaster& Cumulative_uplift);

  /// @brief This updates the cumulative uplift
  /// @param Cumulative_uplift the raster containing the cumulative uplift of this step in m
  ///  passed by reference as this raster is updated in the function
  /// @param Uplift_raster the uplift for the timestep
  /// @param timestep The timestep in years
  /// @date 07/05/2021
  /// @author SMM
  void update_cumulative_uplift(LSDRaster& Cumulative_uplift, LSDRaster& Uplift_raster, float timestep);

  /// @brief Looks at another raster, checks to see if it the same dimensions as the 
  ///  model data, and then replaces any pixel in the model data with nodata
  ///  if it is nodata in the other raster. 
  /// @param OtherRaster the other raster
  /// @author SMM
  /// @date 9/10/2020
  void mask_to_nodata_impose_nodata(LSDRaster& OtherRaster);


  /// @brief create a normal fault horizontally across the raster.
  /// @author FJC
  /// @date 06/07/18
  void normal_fault_part_of_raster(int throw_amt, int throw_type);

  /// @brief simulate instantaneous base level fall
  /// @author FJC
  /// @date 18/07/18
  void base_level_fall(int uplift_amt);

  /// @brief This changes the elevation of a raster
  /// @param elevation_adjust The change in elevation
  /// @author SMM
  /// @date 30/01/2020
  void AdjustElevation(float elevation_change);

  /// @brief This initialises a surface with a hillslope
  /// that is the solution to the nonlinear sediment flux equation.
  /// It overwrites RasterData
  /// @details The parameters D and S_c are stored as data members
  /// Solution from Roering et al., (EPSL, 2007)
  /// @param U the uplift rate
  /// @author SMM
  /// @date 01/07/2014
  void initialise_nonlinear_SS(float U);

  /// @brief This resizes the LSDRasterModel, resetting some flags in the process,
  /// as well as setting many of the Array2D data members to be empty arrays
  /// The raster data in the end is a random surface (determined by the noise
  /// data member)
  /// @param new_rows the new number of rows
  /// @param new_cols the new number of columns
  /// @author SMM
  /// @date 30/06/2014
  void resize_and_reset( int new_rows, int new_cols );

  /// @brief This resizes the LSDRasterModel, resetting some flags in the process,
  /// as well as setting many of the Array2D data members to be empty arrays
  /// The raster data in the end is a random surface (determined by the noise
  /// data member)
  /// This overloaded version also changes the data resolution
  /// @param new_rows the new number of rows
  /// @param new_cols the new number of columns
  /// @param new_resolution the new data resolution
  /// @author SMM
  /// @date 30/06/2014
  void resize_and_reset( int new_rows, int new_cols, float new_resolution );

  /// @brief this ads a pathname to the default names
  /// @param the name of the path
  /// @author SMM
  /// @date 18/06/2014
  void add_path_to_names( string pathname);


  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  // TOOLS FOR TRANSIENT RUNS
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief This takes some "phases" with rates and calculates the rate
  ///  for a given phase of incision
  /// @param phase_start A vector or starting times for each phase
  /// @param phase_rates A vector of the base level drop rates
  /// @return the rate for a given time (the time is withing the data member current_time)
  /// @author SMM
  /// @date 08/02/2021
  float calculate_bl_drop_rate( vector<float> phase_start, vector<float> phase_rates);

  /// @brief This takes some "phases" with rates and calculates the rate
  ///  for a given phase of incision
  /// @param phase_start A vector or starting times for each phase
  /// @param phase_rates A vector of vectors with drop rates at specific pixels
  /// @return a vector of rates for a given time (the time is within the data member current_time)
  /// @author SMM
  /// @date 23/04/2021
  vector<float> calculate_bl_drop_rate( vector<float> phase_start, vector< vector<float> > phase_rates);

  /// @brief This takes some "phases" with elevations calculates the rate
  ///  for a given phase of incision
  /// @param phase_start A vector or starting times for each phase
  /// @param phase_elevations A vector of vectors with elevations of the base level
  /// @return a rate for a given time (the time is within the data member current_time)
  /// @author SMM
  /// @date 07/08/2021
  float calculate_bl_drop_rate_from_elevations( vector<float> phase_start, vector<float> phase_elevations);

  /// @brief This takes some "phases" with elevations and calculates the rates between those phases
  /// @param phase_start A vector or starting times for each phase
  /// @param phase_rates A vector of vectors with elevations at specific pixels
  /// @return a vector of vectors with rates that correspond to the phases 
  /// @date 23/04/2021
  vector< vector<float> > calculate_phase_rates_from_elevations( vector<float> phase_start, vector< vector<float> > phase_elevations);





  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  // TOOLS FOR CHECKING STEADY STATE
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief This function checks to see if the model has achieved steady state.
  /// The nature of steady state checked is set by the cycle_steady_check flag
  /// If this is false, it checks if there is simple steady state
  /// (the surface elevations do not change in time)
  /// If the cycle_steady_check is true, it check if steady state has been
  /// achieved from one cycle to another
  /// @return does not return anything, but instead changes the steady_state flag
  /// @author JAJ     commented SMM
  /// @date 01/01/2014 commented SMM 27/06/2014
  void check_steady_state( void );

  /// @brief This function checks to see if the model should record results
  /// If initial steady state has not been reached, recording is set to
  /// false: that is, the  model does not record information on the build up
  /// to steady state.
  /// @author JAJ
  /// @date 01/01/2014
  void check_recording( void );

  /// @brief This checks on the ending condition of the model run
  /// endTime_mode:
  ///  1 == The end time is just some fixed time after initial steady state
  ///  2 == The end time is after a fixed number of cycles
  ///  3 == The time is after steady state, but waits for a fixed number of cycles
  ///       before ending
  /// @return returns a boolean that is true if the end time has been reached
  /// and false if end time has not been reached
  /// @author JAJ, comments SMM
  /// @date 01/01/2014 comments SMM 27/06/2014
  bool check_end_condition( void );

  /// @brief This function checks to see if this is a periodic run. If it is,
  /// it sets the times to align with the period
  /// Note this only realy comes into play if period mode == 2 or 4
  /// period_mode means
  /// 1 (default) one periodicity used without
  /// 2 Two periodicities that switch at a given interval
  /// 3 Two periodicities used as a compound sin wave
  /// 4 Same as three, but weightings switch at a given interval (as in 2)
  /// @author JAJ
  /// @date 01/012014
  void check_periodicity_switch( void );

  /// @brief If the periodic model cycles over 100 times this returns true
  /// @author JAJ
  bool check_if_hung( void );

  /// @brief Reset model - reset erosion values to 0 after a complete model run
  /// @author JAJ
  /// @date 01/01/2014
  void reset_model( void );

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // @@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@
  // Deal with the BOUNDARY CONDITIONS
  // @@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief This first function is used as a simple way to implement boundary conditions,
  /// particularly no flux and periodic boundary conditions.
  /// The buffered surface has NRows+2 rows and NCols+2 columns.
  /// The integer b_type sets the type of boundary conditions, but currently there
  /// is only one implementation: no flux across N and S; periodic for E and W.
  /// @param b_type at the moment this is irrelevant since this just switches to default
  /// @return creates a buffered LSDRasterModel
  LSDRasterModel create_buffered_surf(int b_type);

  /// @brief This second version has periodic boundaries at E and W boundaries, and
  /// Neumann boundary conditions (prescribed elevations) at the N and S
  /// boundaries.
  /// @param South_boundary_elevation the elevation at the southern boundary
  /// @param North_boundary_elevation the elevation at the southern boundary
  LSDRasterModel create_buffered_surf(float South_boundary_elevation,float North_boundary_elevation);

  /// @brief Check whether current node is a base level node
  /// @param row
  /// @param column
  /// @return true or false
  /// @author JAJ
  /// @date 01/01/2014
  bool is_base_level(int i, int j);

  /// @brief not sure what this does yet (SMM)
  void interpret_boundary(short &dimension, bool &periodic, int &size);

  /// @brief Gets the maxium elevation along a boundary
  /// @param boundary_number 0 == row 0   1 == col 0
  /// @return the maximum elevation along the boundaty
  float find_max_boundary( int boundary_number );


  /// @brief This fixes a channel, derived from source points data
  ///  onto the model DEM
  /// @param source_points_data an LSDSpatialCSVReader object. It needs lat and long and elevation columns
  /// @param column_name the name of the elevation column
  /// @author SMM
  /// @date 04/03/2020
  void impose_channels(LSDSpatialCSVReader& source_points_data, string column_name);

  /// @brief This fixes a channel, derived from source points data
  ///  onto the model DEM
  /// @param source_points_data an LSDSpatialCSVReader object. It needs lat and long and elevation columns
  /// @param column_name the name of the elevation column
  /// @param bl_code a code for baselevel nodes. Usually there is a baselevel code that indicates the node is to be ignored
  /// @author SMM
  /// @date 04/03/2020
  void impose_channels(LSDSpatialCSVReader& source_points_data, string column_name, vector<int> bl_code);


  /// @brief This fixes a channel, derived from source points data
  ///  onto the model DEM. It also lifts the raster so no pixels in the raster are lower than
  /// @param source_points_data an LSDSpatialCSVReader object. It needs lat and long and elevation columns
  /// @param tcolumn_name he name of the elevation column
  /// @author SMM
  /// @date 09//2020
  void impose_channels_and_lift_raster(LSDSpatialCSVReader& source_points_data, string column_name);

  /// @brief Takes a model step and gets an LSDSpatialCSVReader object for later use
  /// @param contributing_pixels for the channel network
  /// @return source_points_data an LSDSpatialCSVReader object. It needs lat and long and elevation columns
  /// @author SMM
  /// @date 04/03/2020
  LSDSpatialCSVReader  get_channels_for_burning(int contributing_pixels, string temp_channel_path, string temp_channel_fname);

  /// @brief Caps elevations using the initial raster
  ///  WARNING no testing if the raster is the correct shape!
  /// @param InitialRaster The initial raster above which the new surface cannot rise.
  /// @author SMM
  /// @date 31/08/2020
  void cap_elevations(LSDRaster& InitialRaster);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // @@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@
  // Calculate erosion rates
  // @@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief Simple function that creates an array with the erosion rates for a given
  /// timestep. It doesn't do anything with NoData cells, and for cells with data
  /// it calls the get_erosion_at_cell member function
  /// updated to catch instances when zeta_old has not been calculated
  /// @return an array with the erosion rates
  /// @author JAJ  updated SMM
  /// @date 01/01/2014  updated 01/07/2014
  Array2D<float> calculate_erosion_rates( void );

  /// @brief Simple function that creates an array with the erosion rates for a given
  /// timestep. It doesn't do anything with NoData cells, 
  /// and uses an uplift raster
  /// @param Uplift_raster a raster of uplift values in m/yr
  /// @return an array with the erosion rates
  /// @author SMM
  /// @date 6/05/2021
  Array2D<float> calculate_erosion_rates( LSDRaster& Uplift_raster );

  /// @brief This calculates the erosion rate for individual cells.
  /// Currently it assumes that the zeta_old data member is from the previous
  /// timestep
  /// @param row the row of the cell
  /// @param col the column of the cell
  /// @author JAJ
  /// @date 01/01/2014
  float get_erosion_at_cell(int row, int col);

  /// @brief this calcualtes the total erosion over a timester
  /// @return the erosion rate calculated over the last timestep
  /// @author SMM
  /// @date 01/08/2014
  float get_total_erosion_rate_over_timestep();

  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  /// CREATE PRECIPITION FLUX ARRAY
  /// Produces precipitation array from provided precipitation rate.
  ///---------------------------------------------------------------------------
  Array2D<float> precip_array_from_precip_rate(float precip_rate);

  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// TOPOGRAPHIC DERIVATIVES
  /// Specifically, this function gets the topographic slopes, as required for
  /// the sediment flux calculations.  The slopes are stored as two matrices, one
  /// that stores slopes between rows, the other which for slopes between
  /// columns.  Note that this is a finite volume model that utilises cubic model
  /// voxels. Sediment fluxes are only permitted through the faces.
  ///
  /// For slopes between columns, the entry at S[row][col] refers to the slope
  /// between zeta at node [row][col] and at node [row][col+1].  Likewise for the
  /// slopes between rows.  In short, the center points of the slopes are offset
  /// by 1/2 a node spacing in the positive direction.
  /// Note that there are NCols +1 and NRows +1 columns and rows respectively
  ///----------------------------------------------------------------------------
  void get_slopes(Array2D<float>& SlopesBetweenRows, Array2D<float>& SlopesBetweenCols);
  ///----------------------------------------------------------------------------
  /// get_topographic_divergence
  /// gets the topographic divergence at each point in the model domain.  Use
  /// buffered topography
  ///----------------------------------------------------------------------------
  Array2D<float> get_topographic_divergence();

  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// HYDROLOGICAL TOOLS
  ///----------------------------------------------------------------------------
  /// calculate_channel_width_wolman
  /// This function calculates channel width using the wolman method.
  /// NOTE: typically Q_w will be in m^3/s.
  /// EXAMPLE: in Salmon River, Idaho (Emmett, 1975 cited in Knighton 1988):
  ///          k_w = 2.77 and b = 0.56. b is often assumed to be 0.5
  ///----------------------------------------------------------------------------
  float calculate_channel_width_wolman(float Q_w, float k_w, float b);
  ///----------------------------------------------------------------------------
  /// array_channel_width_wolman
  /// this function calcualtes channel width in a stand alone module so the widths
  /// can be tested
  ///----------------------------------------------------------------------------
  Array2D<float> array_channel_width_wolman(Array2D<float>& Q_w, float& k_w, float& b);

  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// EROSION RATES/SEDIMENT FLUXES
  ///
  ///----------------------------------------------------------------------------
  /// this caluclates the fluvial erosion rate at each point
  ///----------------------------------------------------------------------------
  Array2D<float> calculate_fluvial_erosion_rate(Array2D<float> ChannelWidth, Array2D<float> Q_w,
    Array2D<float> TopoDivergence, float K, float n, float m, float eros_thresh);

  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// IMPLICIT MODEL COMPONENTS
  ///------------------------------------------------------------------------------
  /// Implicit schemes for combination of hillslope sediment transport using
  /// non-linear hillslope transport law, and fluvial erosion.  This is essentially
  /// the implicit implementation of MuddPILE, but has been modified so that now
  /// fluvial erosion is undertaken using FASTSCAPE (Braun and Willet, 2013), which
  /// greatly increases computational efficiency.
  ///------------------------------------------------------------------------------
  /// calculate_k_values_for_assembly_matrix/mtl_initiate_assembler_matrix
  /// this function creates vectors of integers that refer to the k values, that is
  /// the index into the vectorized matrix of zeta values, that is used in the assembly matrix
  /// the number of elements in the k vectors is N_rows*N_cols
  ///------------------------------------------------------------------------------
  void calculate_k_values_for_assembly_matrix(int NRows, int NCols, vector<int>& k_value_i_j,
                      vector<int>& k_value_ip1_j,  vector<int>& k_value_im1_j, vector<int>& k_value_i_jp1,
                      vector<int>& k_value_i_jm1);

  /// mtl_initiate_assembler_matrix
  void mtl_initiate_assembler_matrix(int& problem_dimension,
           float& inv_dx_S_c_squared, float& inv_dy_S_c_squared, float& dx_front_term,
               float& dy_front_term, vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,
           vector<int>& vec_k_value_im1_j, vector<int>& vec_k_value_i_jp1,
            vector<int>& vec_k_value_i_jm1);
  //------------------------------------------------------------------------------
  /// mtl_assemble_matrix
  /// this function assembles the solution matrix for nonlinear creep transport
  //------------------------------------------------------------------------------
  void mtl_assemble_matrix(Array2D<float>& zeta_last_iter, Array2D<float>& zeta_last_timestep,
    Array2D<float>& zeta_this_iter, Array2D<float>& uplift_rate,
                Array2D<float>& fluvial_erosion_rate,
             mtl::compressed2D<float>& mtl_Assembly_matrix, mtl::dense_vector<float>& mtl_b_vector,
     float dt, int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
             float dx_front_term, float dy_front_term,
             float South_boundary_elevation, float North_boundary_elevation,
             vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,vector<int>& vec_k_value_im1_j,
     vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1);
  ///------------------------------------------------------------------------------
  /// mtl_solve_assembler_matrix
  /// this function assembles the solution matrix
  ///------------------------------------------------------------------------------
  void mtl_solve_assembler_matrix(Array2D<float>& zeta_last_iter, Array2D<float>& zeta_last_timestep,
    Array2D<float>& zeta_this_iter, Array2D<float>& uplift_rate,
                Array2D<float>& fluvial_erosion_rate,
    float dt, int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
             float dx_front_term, float dy_front_term,
             vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j, vector<int>& vec_k_value_im1_j,
             vector<int>& vec_k_value_i_jp1, std::vector<int>& vec_k_value_i_jm1,
             float South_boundary_elevation, float North_boundary_elevation);
  ///------------------------------------------------------------------------------
  /// nonlinear_creep_timestep
  /// do a creep timestep.  This function houses the above two functions to
  /// undertake model timestep using implicit implementation of the nonlinear
  /// transport law.
  /// NOTE you need to run mtl_initiate_assembler_matrix before you run this function
  ///------------------------------------------------------------------------------
  void nonlinear_creep_timestep(Array2D<float>& fluvial_erosion_rate, float iteration_tolerance,
        int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
        float dx_front_term, float dy_front_term, vector<int>& vec_k_value_i_j,
        vector<int>& vec_k_value_ip1_j, vector<int>& vec_k_value_im1_j,
        vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1,
        float South_boundary_elevation, float North_boundary_elevation);

  /// -----------------------------------------------------------------------------
  /// Soil diffusion method
  /// Container for all the finite volume components
  /// -----------------------------------------------------------------------------
  void soil_diffusion_fv( void );

  /// -----------------------------------------------------------------------------
  /// Finite difference matrix
  /// -----------------------------------------------------------------------------
  mtl::compressed2D<float> generate_fd_matrix( int dimension, int size, bool periodic );
  mtl::dense_vector <float> build_fd_vector( int dimension, int size );
  //void repack_fd_vector(mtl::dense_vector <float> &data_vector, int dimension);

  mtl::compressed2D<float> generate_fv_matrix( int dimension, int size, bool periodic );
  mtl::dense_vector <float> build_fv_vector( int dimension, int size );
  void repack_vector(mtl::dense_vector <float> &data_vector, int dimension);

  /// -----------------------------------------------------------------------------
  /// Soil diffusion using linear flux model
  /// Solved using finite difference
  /// -----------------------------------------------------------------------------
  void soil_diffusion_fd_linear( void );

  void soil_diffusion_fv_nonlinear( void );

  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// RUN MODEL
  ///------------------------------------------------------------------------------
  /// A series of wrapper functions that implement the numerical model
  ///------------------------------------------------------------------------------
  /// implicit_hillslope_and_fluvial
  /// This function sets up a landscape evolution model run incorporating fluvial
  /// erosion and hillslope erosion via non-linear creep.  It calls the implicit
  /// implementation and returns the topography after the final timestep.
  /// The user should provide the parameter file which sets out the details of the
  /// model run.
  ///------------------------------------------------------------------------------
  LSDRasterModel run_model_implicit_hillslope_and_fluvial(string param_file);

  /// @brief This wrapper just calls the run_components method. Parameters used are those
  /// stored as data members
  /// @author JAJ
  /// @date 01/01/2014
  void run_model( void );

  /// @brief This loads a steady state raster from the data members so that
  /// a model can be run repeatedly from the same steady state condition
  /// @author JAJ
  /// @date 01/01/2014
  void run_model_from_steady_state( void );

  /// @brief This wrapper just calls the run_components method. Parameters used are those
  /// stored as data members. This one actually calls the erosion laws
  /// @author JAJ
  /// @date 01/01/2014
  void run_components( void );

  /// @brief This is a wrapper similar to run_components but sends the
  /// fluvial and uplfit fields to the nonlinear solver
  /// @author SMM
  /// @date 07/07/2014
  void run_components_combined( void );

  /// @brief This is a wrapper similar to run_components but sends the
  /// fluvial and uplfit fields to the nonlinear solver.
  /// @detail Variable U and K rasters can be used.
  /// @param URaster A raster of uplift rates
  /// @param KRaster A raster of K values
  /// @param use_adaptive_timestep If true, an adaptive timestep is used
  /// @author SMM
  /// @date 03/09/2017
  void run_components_combined( LSDRaster& URaster, LSDRaster& KRaster, bool use_adaptive_timestep );

  /// @brief This is a wrapper similar to run_components but sends the
  /// fluvial and uplfit fields to the nonlinear solver.
  /// This one allows you to put in a base level
  /// And include transience
  /// @detail Variable U and K rasters can be used.
  /// @param URaster A raster of uplift rates
  /// @param KRaster A raster of K values
  /// @param use_adaptive_timestep If true, an adaptive timestep is used
  /// @param source_points_data some points to serve as base level nodes
  /// @param phase_time a vector of times for different base level fall rates
  /// @param phase_rates a vector of base level fall rates
  /// @param column_name the name of the elevation column
  /// @author SMM
  /// @date 08/02/2021
  void run_components_combined_imposed_baselevel( LSDRaster& URaster, LSDRaster& KRaster, 
                                          bool use_adaptive_timestep, LSDSpatialCSVReader& source_points_data,
                                          vector<float> phase_time, vector<float> phase_rates, string column_name);


  /// @brief This is a wrapper similar to run_components but sends the
  /// fluvial and uplfit fields to the nonlinear solver.
  /// This one allows you to put in a base level
  /// And include transience. This one has transience at mulitiple pixels
  /// @detail Variable U and K rasters can be used.
  /// @param URaster A raster of uplift rates
  /// @param KRaster A raster of K values
  /// @param use_adaptive_timestep If true, an adaptive timestep is used
  /// @param source_points_data some points to serve as base level nodes
  /// @param phase_time a vector of times for different base level fall rates
  /// @param phase_eleations a vec vec of base level elevations from which rates are calculated. 
    /// @param minimum_slope the minimum slope for elevation change between pixels. This replaces overexcavated nodes
  /// @param column_name the name of the elevation column
  /// @author SMM
  /// @date 23/04/2021
  void run_components_combined_imposed_baselevel( LSDRaster& URaster, LSDRaster& KRaster, 
                                          bool use_adaptive_timestep, LSDSpatialCSVReader& source_points_data,
                                          vector<float> phase_time, vector< vector<float> > phase_elevations, 
                                          float minimum_slope, string column_name);


  /// @brief This is a wrapper similar to run_components but sends the
  /// fluvial and uplfit fields to the nonlinear solver.
  /// This one allows you to put in a base level
  /// And include transience. This one has transience at mulitiple pixels
  /// It also has the lithocube. 
  /// It is a beast, really
  /// @detail Variable U and K rasters can be used.
  /// @param URaster A raster of uplift rates
  /// @param use_adaptive_timestep If true, an adaptive timestep is used
  /// @param source_points_data some points to serve as base level nodes
  /// @param phase_time a vector of times for different base level fall rates
  /// @param phase_eleations a vec vec of base level elevations from which rates are calculated. 
  /// @oaram cumulative_uplift a raster of cumulative uplift for exhumation calculations
  /// @param LSDLC a lithocube object
  /// @param forbidden_lithocodes a list of lithocodes that will be set to the default value
  /// @param print_lithocode_raster a bool when true print the lithocode raster at the printing timesteps
  /// @param use_hillslope_hybrid if true, turns on the hybrid model that does linear diffusion but uses a critical
  ///  slope if too steep
  /// @param threshold_contributing_pixels Number of pixels to designate as hillslopes for the hillslope model
  /// @param minimum_slope the minimum slope for elevation change between pixels. This replaces overexcavated nodes
  /// @param print_exhumation_and_cumulative_uplift if true prints the exhumation surface and cumulative uplift
  /// @param column_name the name of the elevation column
  /// @author SMM
  /// @date 06/05/2021
  void run_components_combined_imposed_baselevel( LSDRaster& URaster, 
                                          bool use_adaptive_timestep, LSDSpatialCSVReader& source_points_data,
                                          vector<float> phase_time, vector< vector<float> > phase_elevations,
                                          LSDRaster& cumulative_uplift, 
                                          LSDLithoCube& LSDLC, list<int> forbidden_lithocodes, 
                                          bool print_lithocode_raster,
                                          bool use_hillslope_hybrid,
                                          int threshold_contributing_pixels, 
                                          float minimum_slope,
                                          bool print_exhumation_and_cumulative_uplift,
                                          string column_name);


  /// @brief This is a wrapper similar to run_components but sends the
  /// fluvial and uplfit fields to the nonlinear solver.
  /// This one allows you to put in a base level
  /// And include transience. This one has transience at mulitiple pixels
  /// It also has the lithocube. 
  /// It is a beast, really
  /// @detail Variable U and K rasters can be used.
  /// @param URaster A raster of uplift rates
  /// @param use_adaptive_timestep If true, an adaptive timestep is used
  /// @param source_points_data some points to serve as base level nodes
  /// @param phase_time a vector of times for different base level fall rates
  /// @param outlet_elevations is a vector of elevations at the outlet 
  /// @oaram cumulative_uplift a raster of cumulative uplift for exhumation calculations
  /// @param LSDLC a lithocube object
  /// @param forbidden_lithocodes a list of lithocodes that will be set to the default value
  /// @param print_lithocode_raster a bool when true print the lithocode raster at the printing timesteps
  /// @param use_hillslope_hybrid if true, turns on the hybrid model that does linear diffusion but uses a critical
  ///  slope if too steep
  /// @param threshold_contributing_pixels Number of pixels to designate as hillslopes for the hillslope model
  /// @param minimum_slope the minimum slope for elevation change between pixels. This replaces overexcavated nodes
  /// @param print_exhumation_and_cumulative_uplift if true prints the exhumation surface and cumulative uplift
  /// @param column_name the name of the elevation column
  /// @author SMM
  /// @date 07/08/2021
  void run_components_combined_imposed_baselevel( LSDRaster& URaster, 
                                          bool use_adaptive_timestep, LSDSpatialCSVReader& source_points_data,
                                          vector<float> phase_time, vector< float > outlet_elevations,
                                          LSDRaster& cumulative_uplift, 
                                          LSDLithoCube& LSDLC, list<int> forbidden_lithocodes, 
                                          bool print_lithocode_raster,
                                          bool use_hillslope_hybrid,
                                          int threshold_contributing_pixels, 
                                          float minimum_slope,
                                          bool print_exhumation_and_cumulative_uplift,
                                          string column_name);





  /// @brief This is a wrapper that runs the model but includes CRN columns
  /// fluvial and uplfit fields to the nonlinear solver
  /// @param CRNColumns the vector of particle columns
  /// @param eroded_cells this gets replaced, it is the eroded particles
  /// @param startType the starting type of the particle (doesn't really play a role)
  /// @param startDepth the starting depth (in m) of the particles
  /// @param particle_spacing vertical distance between particles in the column
  /// @param CRNParam the cosmogenic parameter values
  /// @author SMM
  /// @date 25/07/2014
  void run_components_combined_cell_tracker( vector<LSDParticleColumn>& CRNColumns,
                      vector<LSDParticleColumn>& eroded_cells,
                      int startType, double startDepth, double particle_spacing,
                      LSDCRNParameters& CRNParam);

  /// @brief This initiates a vector of CRN columns that sit under the model
  /// @author SMM
  /// @param column_spacing and integer telling how many nodes between particle columns
  /// @param CRNcol_rows this is an integer vector that is replaced in this function
  /// each element indexes the row of the vector of columns
  /// @param CRNcol_cols this is an integer vector that is replaced in this function
  /// each element indexes the col of the vector of columns
  /// @param rho_r the density of the rock in kg/m^3
  /// @param this_U the uplift rate in m/yr that the particles CRN concentrations
  /// will be equilibrated to. Note it is only via nucleonic production
  /// @param startType the starting type of the particle (doesn't really play a role)
  /// @param startDepth the starting depth (in m) of the particles
  /// @param particle_spacing vertical distance between particles in the column
  /// @param CRNParam the cosmogenic parameter values
  /// @return a vector of particle columns
  /// @author SMM
  /// @date 31/07/2014
  vector<LSDParticleColumn> initiate_steady_CRN_columns(int column_spacing,
                    vector<int>& CRNcol_rows, vector<int>& CRNcol_cols,
                    double rho_r, double this_U, int startType, double startDepth,
                    double particle_spacing, LSDCRNParameters& CRNParam);


  /// @brief This method forces the landscape into its steady state profile, by using periodic forcing.
  /// This is much more efficient than using static forcing (as in run model), but doesn't give
  /// a nice animation of an evolving landscape
  /// Swings and roundabouts
  /// @author JAJ
  /// @date 01/01/2014
  void reach_steady_state( void );

  /// @brief This method creates a steady landscape that assumes everywhere obeys the
  ///  stream power law. It is based on equation 4a from Mudd et al 2014 JGR-ES
  /// @detail There is no return but the underlying raster data will reflect the
  ///  analytical steady topography for the uplift rate.
  /// @param U the uplift rate (in m/yr)
  /// @author SMM
  /// @date 10/08/2017
  void fluvial_snap_to_steady_state(float U);

  /// @brief This method creates a steady landscape that assumes everywhere obeys the
  ///  stream power law. It is based on equation 4a from Mudd et al 2014 JGR-ES.
  ///  The function is given a target relief and the K value is adjusted to match
  ///  this target relief at the maximum chi value.
  /// @detail In addition to the return, the underlying raster data will reflect the
  ///  analytical steady topography for the uplift rate.
  /// @param U the uplift rate (in m/yr)
  /// @param desired_relief The desired landscape relief in metres
  /// @return The back calculated K value for the desired relief
  /// @author SMM
  /// @date 10/08/2017
  float fluvial_snap_to_steady_state_tune_K_for_relief(float U, float desired_relief);

  /// @brief This method calcualtes the fluvial K required to generate the
  ///   desired relief at steady state for the farthest upstream chi.
  ///   It is based on equation 4a from Mudd et al 2014 JGR-ES.
  /// @detail This does not update aything in the model, but simply returns the desired K
  /// @param U the uplift rate (in m/yr)
  /// @param desired_relief The desired landscape relief in metres
  /// @return The back calculated K value for the desired relief
  /// @author SMM
  /// @date 29/08/2017
  float fluvial_calculate_K_for_steady_state_relief(float U, float desired_relief);

  /// @brief This smooths the model DEM and returns this smoothed surface as a raster
  /// @param central_pixel_weighting A float that give the weighting of the central pixel. 
  ///   The higher this number, the less smoothing. 2 is probably a good starting value. 
  /// @return The smoothed elevation
  /// @author SMM
  /// @date 16/12/2019
  LSDRaster basic_smooth(float central_pixel_weighting);

  /// @brief Calucaltes the laplacian on the surface. Can be used in the hillslope module
  /// @return The laplacian curvature as a raster
  /// @author SMM
  /// @date 07/05/2021
  LSDRaster basic_curvature();

  /// @brief This checks for rivers (using a drainage area threshold) and then any remaining pixels
  ///  are popped to a critical slope. Creates a river network with striaght slopes in between
  /// @detail Very rudimentary: only uses slopes in the D8 direction so the slopes will 
  ///  not be very accurate if the polyfit slope code is run. It should be considered a maximum relief
  /// @param critical_slope the critical gradient for each D8 connection between pixels
  /// @param contributing_pixel_threshold the threshold in pixles for a channel to form
  /// @author BG, edited by SMM
  /// @date 09/10/2019
  LSDRaster basic_valley_fill_critical_slope(float critical_slope, int contributing_pixel_threshold);

  /// @brief This checks for rivers (using a drainage area threshold) and then any remaining pixels
  ///  are popped to a critical slope. Creates a river network with striaght slopes in between
  /// @brief Very rudimentary: only uses slopes in the D8 direction so the slopes will 
  ///  not be very accurate if the polyfit slope code is run. It should be considered a maximum relief
  /// @param S_c_raster a raster with S_c values. Need to be same dimensions as the model raster
  /// @param contributing_pixel_threshold the threshold in pixles for a channel to form
  /// @author SMM
  /// @date 10/10/2019
  LSDRaster basic_valley_fill_critical_slope(LSDRaster& S_c_raster, int contributing_pixel_threshold);

  /// @brief This method snaps to steady with spatially variable uplift and erodibility fields
  ///  It also allows fixed base level. 
  /// @param K_values a raster of erodiblity
  /// @param U_values a raster of uplift
  /// @param Source_points_data a spatialc csv reader with the appropriate file
  /// @param carve_before_fill if true, run the carving algorithm before the filling algorithm
  /// @param column_name the name of the elevation column
  /// @author SMM
  /// @date 01/10/2019
  void fluvial_snap_to_steady_variable_K_variable_U(LSDRaster& K_values, LSDRaster& U_values, 
                                                    LSDSpatialCSVReader& Source_points_data, 
                                                    bool carve_before_fill, string column_name);


  /// @brief This is a hybrid model that calculates hillslope diffusion as well as a critical slope
  ///  so the diffused hillslopes cannot exceed a critical slope 
  /// @detail Must be used after the fluvial step
  /// @param Sc_values a raster of critical slopes
  /// @param U_values a raster of uplift
  /// @param threshold_pixels The number of contributing pixels must be less than this to snap the hillslope
  /// @param Source_points_data a list of baselevel nodes to ignore
  /// @author SMM
  /// @date 11/05/2021
  void hillslope_hybrid_module(LSDRaster& Sc_values, LSDRaster& Uplift, 
                               int threshold_pixels,LSDSpatialCSVReader& Source_points_data);


  /// @brief A hillslope snapping routine based on the stack. 
  /// @param Sc_values a raster of critical slopes
  /// @param carve_before_fill if true, run the carving algorithm before the filling algorithm
  /// @param threshold_pixels The number of contributing pixels must be less than this to snap the hillslope
  /// @author SMM
  /// @date 06/05/2021
  void hillslope_snap_to_steady_variable_Sc(LSDRaster& Sc_values, bool carve_before_fill, int threshold_pixels);

  /// @brief This method instantaneously tilts the landscape by a certain angle.
  /// @param angle the tilt angle in degrees
  /// @param tilt_boundary. Tilt can be from the N, E, S, or W boundary. Must be"N", "E",
  /// "S", or "W".
  /// @return Doesn't return anlything but updates the raster elevations
  /// @author FJC
  /// @date 11/03/19
  void instantaneous_tilt(float angle, string tilt_boundary);


  /// @brief Fastscape, implicit finite difference solver for stream power equations
  /// O(n)
  /// Method takes its paramaters from the model data members
  /// and solves the stream power equation at a future timestep in linear time
  /// @author JAJ, commented SMM
  /// @date 01/01/2014, edit 18/01/2014
  void fluvial_incision( void );

  /// @brief Fastscape, implicit finite difference solver for stream power equations
  /// O(n)
  /// Method takes its paramaters from the model data members
  /// and solves the stream power equation at a future timestep in linear time
  /// This version includes the current uplift, so you do not need to call
  /// uplift after this has finished
  /// @author SMM
  /// @date 7/07/2014
  void fluvial_incision_with_uplift( void );

  /// @brief Fastscape, implicit finite difference solver for stream power equations
  /// O(n)
  /// Method takes the K value from a raster fed to it
  /// and solves the stream power equation at a future timestep in linear time
  /// This version includes the current uplift, so you do not need to call
  /// uplift after this has finished
  /// @param K_raster the raster of K values.
  /// @author SMM
  /// @date 01/09/2017
  void fluvial_incision_with_uplift_and_variable_K( LSDRaster& K_raster );

  /// @brief Fastscape, implicit finite difference solver for stream power equations
  /// O(n)
  /// Method takes the K value from a raster fed to it
  /// and also take a raster of the uplift rates
  /// and solves the stream power equation at a future timestep in linear time
  /// This version includes the current uplift, so you do not need to call
  /// uplift after this has finished
  /// @param K_raster the raster of K values.
  /// @param Uplift_rate a raster of uplift rates in m/yr
  /// @author SMM
  /// @date 01/09/2017
  void fluvial_incision_with_variable_uplift_and_variable_K( LSDRaster& Uplift_rate, LSDRaster& K_raster );


  /// @brief Fastscape, implicit finite difference solver for stream power equations
  /// O(n)
  /// Method takes the K value from a raster fed to it
  /// and also take a raster of the uplift rates
  /// and solves the stream power equation at a future timestep in linear time
  /// This version includes the current uplift, so you do not need to call
  /// uplift after this has finished. Uses an adaptive timestep.
  /// @param K_raster the raster of K values.
  /// @param Uplift_rate a raster of uplift rates in m/yr
  /// @author SMM
  /// @date 06/09/2017
  void fluvial_incision_with_variable_uplift_and_variable_K_adaptive_timestep( LSDRaster& Uplift_rate, LSDRaster& K_raster );

  /// @brief This is used to process a transient base level file
  /// @param source_points_data a spatial csv object that has the correct elevation input
  /// @param timing_prefix the prefix of a time column in the csv
  /// @param timing_multiplier how many years are in the each number of timestep (so, say t5 = 5000 year, then the 
  ///   multiplier is 1000)
  /// @param n_time_columns The number of time columns to attempt to read
  /// @param phase_steps How many steps there are between phases. So if you have t5, then t10, then t15 the steps would be 5
  /// @author SMM
  /// @date 22/04/2021
  void process_baselevel( LSDSpatialCSVReader& source_points_data, string timing_prefix, 
                                        float timing_multiplier, int n_time_columns, int phase_steps,
                                        vector<float>& phases_vec, vector< vector<float> >& elevation_vecvec);


  /// @brief This is used to process a transient base level file using only a single elevation value
  /// @param transient_infile_name name of the transient input file including csv
  /// @param phases Replaced in code these are the times
  /// @param outlet_elevations these are the elevations at the fixed times of the outlet
  /// @author SMM
  /// @date 07/08/2021
  void process_transient_file(string transient_infile_name, vector<float>& phases,vector<float>& outlet_elevations);


  /// @brief This is a little helper to tag specific nodes with a baselelvel switch
  /// @detail creates a boolean column where 0 is an area node and 1 is a node that lowers at a fixed rate. 
  /// @param source_points_data a spatial csv object that has the correct elevation input
  /// @param column_name The name of the column to test against a criteria. 
  /// @author SMM
  /// @date 26/04/2021
  void baselevel_run_area_switch( LSDSpatialCSVReader& source_points_data, string column_name);


  /// @brief Fastscape, implicit finite difference solver for stream power equations
  /// O(n)
  /// Method takes the K value from a raster fed to it
  /// and also take a raster of the uplift rates
  /// and solves the stream power equation at a future timestep in linear time
  /// This version includes the current uplift, so you do not need to call
  /// uplift after this has finished. Uses an adaptive timestep.
  /// This version allows you to impose a base level
  /// @param K_raster the raster of K values.
  /// @param Uplift_rate a raster of uplift rates in m/yr
  /// @param source_points_data a spatial csv object that has the correct elevation input
  /// @param the rate the base level is falling. Need a rate rather than a fixed elevation because
  ///  of the adaptive timestep
  /// @param let_timestep_increase a boolean that when true allows the timestep to increase. Make false
  ///   when you need to ensure the timestep lands on a specific time. 
  /// @param column_name the name of the elevation column
  /// @author SMM
  /// @date 08/02/2021
  void fluvial_incision_with_variable_uplift_and_variable_K_adaptive_timestep_impose_baselevel( LSDRaster& Urate_raster, 
                                 LSDRaster& K_raster, LSDSpatialCSVReader& source_points_data,
                                 float bl_fall_rate, bool let_timestep_increase, string column_name);

  /// @brief Fastscape, implicit finite difference solver for stream power equations
  /// O(n)
  /// Method takes the K value from a raster fed to it
  /// and also take a raster of the uplift rates
  /// and solves the stream power equation at a future timestep in linear time
  /// This version includes the current uplift, so you do not need to call
  /// uplift after this has finished. Uses an adaptive timestep.
  /// This version allows you to impose a base level
  /// @param K_raster the raster of K values.
  /// @param Uplift_rate a raster of uplift rates in m/yr
  /// @param source_points_data a spatial csv object that has the correct elevation input
  /// @param the rate the base level is falling. Need a rate rather than a fixed elevation because
  ///  of the adaptive timestep. This has baselevel rate at a vector of baselevel nodes. 
  /// @param let_timestep_increase a boolean that when true allows the timestep to increase. Make false
  ///   when you need to ensure the timestep lands on a specific time. 
  /// @param column_name the name of the elevation columns
  /// @author SMM
  /// @date 08/02/2021
  void fluvial_incision_with_variable_uplift_and_variable_K_adaptive_timestep_impose_baselevel( LSDRaster& Urate_raster, 
                                 LSDRaster& K_raster, LSDSpatialCSVReader& source_points_data, 
                                 vector<float> bl_fall_rate,
                                 bool let_timestep_increase,
                                 bool use_adaptive_timestep,
                                 float minimum_slope,
                                 string column_name);                               

  /// @brief This function is more or less identical to fluvial_incision above, but it
  /// Returns a raster with the erosion rate and takes arguments rather
  /// than reading from data members
  /// @param timestep the time spacing
  /// @param K fluvial erosivity
  /// @param m area exponent
  /// @param n slope exponent
  /// @param boundary a vector of strings cotaining model boundary conditions
  /// @return A raster containing the erosion rate from fluvial processes
  /// @author JAJ
  /// @date 01/01/2014
  LSDRaster fluvial_erosion_rate(float timestep, float K, float m, float n, vector <string> boundary);

  /// @brief This function is more or less identical to fluvial_incision above, but it
  /// Returns an array and takes arguments rather reads from data members
  /// @return A 2d float containing the erosion rate from fluvial processes
  /// @author SMM
  /// @date 07/07/2014
  Array2D<float> fluvial_erosion_rate( void );

  /// @brief This assumes that all sediment transported from rivers into
  /// channels is removed. It checks the raster to see where the channels are
  /// which at this point is determined by a threshold drainage area,
  /// and then removes all the sediment to those pixels
  /// @author JAJ
  /// @date 01/01/2014
  void wash_out( void );

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  // TOOLS FOR UPLIFT
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief Creates uplift field from a set of templates
  /// @param mode specifies a mode of uplift:
  ///   (0) - block uplift
  ///  1   - tilt block
  ///  2   - gaussian
  ///  3   - quadratic
  /// @param second argument is the maximum uplift
  /// @return returns the uplift field that is the same dimensions as the original
  /// raster
  /// @author JAJ
  /// @date 01/01/2014
  Array2D <float> generate_uplift_field( int mode, float max_uplift);

  /// @brief Creates uplift field from a set of templates, parameters are taken from
  /// data members
  /// @return returns the uplift field that is the same dimensions as the original
  /// raster
  /// @author SMM
  /// @date 07/07/2014
  Array2D <float> generate_uplift_field( void );


  /// @brief Gets the uplift value at a given cell
  /// this method is implemented as a memory saving measure, rather than
  /// storing the uplift field in memory
  /// Some methods still implemented still use this uplift field
  /// It's advisable this is changed, otherwise the size of rasters that
  /// can be modelled will be severely reduced
  /// @details    This uses the uplift_mode to determine how uplift is calculated
  ///   (0) - block uplift
  ///  1   - tilt block
  ///  2   - gaussian
  ///  3   - quadratic
  /// 4   - periodic
  /// @param row
  /// @param column
  /// @return the uplift (as a distance rather than rate, uses data member timestep)
  /// @author JAJ  commented SMM
  /// @author 01/01/2014   commented 26/06/2014
  float get_uplift_at_cell(int i, int j);

  /// @brief Gets the uplift rate at a given cell
  /// this method is implemented as a memory saving measure, rather than
  /// storing the uplift field in memory
  /// Some methods still implemented still use this uplift field
  /// It's advisable this is changed, otherwise the size of rasters that
  /// can be modelled will be severely reduced
  /// @details    This uses the uplift_mode to determine how uplift is calculated
  ///   (0) - block uplift
  ///  1   - tilt block
  ///  2   - gaussian
  ///  3   - quadratic
  /// 4   - periodic
  /// @param row
  /// @param column
  /// @return the uplift rate
  /// @author SMM
  /// @date 07/07/2014   commented 26/06/2014
  float get_uplift_rate_at_cell(int i, int j);

  /// @brief This checks to see if the uplift field is consistent with
  ///  the raster dimensions. If not it corrects the dimensions of the uplift field
  /// @author SMM
  /// @date 23/08/2017
  void check_and_correct_uplift_field();


  /// @brief this calcualtes the average uplfit rate over the entire model domain,
  /// excluding the N and S boundaries
  /// @return the average uplift rate in m/yr
  /// @author SMM
  /// @date 01/08/2014
  float get_average_upflit_rate_last_timestep();


  /// @brief Apply uplift field to the raster.  Overloaded function so that the first
  /// simply considers uniform uplift, the second allows user to use a prescribed
  /// uplift fields of greater complexity, for example taking account of fault
  /// geometry.
  /// @details WARNING the returned LSDRasterModel only contains a very small
  /// subset of the data members of the original LSDRasterModel. Implementation
  /// NOT RECOMMENDED!
  /// @param UpliftRate the rate of uplift at that timestep
  /// @param dt the timestep
  /// @param Returns an LSDRasterModel of uplift (SMM: why not just update
  /// the raster directly??)
  /// @author JAJ, comments SMM
  /// @date 01/01/2014  SMM comments 26/06/2014
  LSDRasterModel uplift_surface(float UpliftRate, float dt);

  /// @brief Uplift surface using specified uplift field
  /// uplift field should be specified as an array with the same dimensions as the
  /// elevation raster, permitting non-uniform uplift fields to be applied in the
  /// model.
  /// @details WARNING the returned LSDRasterModel only contains a very small
  /// subset of the data members of the original LSDRasterModel. Implementation
  /// NOT RECOMMENDED!
  /// @param UpliftRate a 2D float array of the uplift. Can be made using the
  /// member function generate_uplift_field
  /// @param dt the timestep
  /// @param an LSDRasterModel object (SMM again, why not update the underlying
  /// data member of surface elevation?)
  /// @author JAJ, comments SMM
  /// @date 01/01/2014  SMM comments 26/06/2014
  LSDRasterModel uplift_surface(Array2D<float> UpliftRate, float dt);

  /// @brief Intrinsic method of uplifting the Raster
  /// Uplift field attribute is incremented onto RasterData itself
  /// There are no parameters, but rather it simply passes upflift to the
  /// get upflift at cell function
  /// @details Uplift is calculated based on data_members max_uplift, timestep and
  /// uplift mode.
  /// This uses the uplift_mode to determine how uplift is calculated
  ///   (0) - block uplift
  ///  1   - tilt block
  ///  2   - gaussian
  ///  3   - quadratic
  /// 4   - periodic
  /// @author JAJ commented SMM 26/06/2014
  /// @date 01/01/2014, commented 26/06/2014
  void uplift_surface( void );

  /// @brief  This just returns the max_uplift data member
  /// NOTE; while this is currently a very trivial, and arguably unecessary method, it should be used
  /// and developed if someone wants to integrate some sort of changing uplift field
  /// @return the data member holding the maximum uplift
  /// @author JAJ
  /// @date 01/01/2014
  float get_max_uplift( void );

  /// @brief This function sets the uplift_field data member as bolck uplift
  /// with a rate of uplift_rate
  /// @param uplift_rate a float of uplift rate, the entire block will uplift
  /// at this rate
  /// @author SMM
  /// @date 03/07/2014
  void set_uplift_field_to_block_uplift(float uplift_rate);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  // TOOLS FOR ISOSTACY
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// --------------------------------------------------------------------
  /// Runs flexural isostatic calculations
  /// Uses fourier filtering method
  /// Pelletier (2008)
  /// --------------------------------------------------------------------
  LSDRasterModel run_isostatic_correction( void );

  /// -------------------------------------------------------------------
  /// Correct for isostasy using Airy model
  /// -------------------------------------------------------------------
  void Airy_isostasy( void );

  /// -------------------------------------------------------------------
  /// Correct for isostasy using flexural model
  /// -------------------------------------------------------------------
  void flexural_isostasy( float alpha );
  void flexural_isostasy_alt( void );

  void write_root(string name, string ext);

  /// -------------------------------------------------------------------
  /// Calculates depth of topographic root using FFT methods inherited from LSDRasterSpectral
  /// -------------------------------------------------------------------
  Array2D <float> calculate_root( void );
  Array2D <float> calculate_airy( void );

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  // Setter methods
  // @~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief this sets the boundary conditions
  void set_boundary_conditions(vector <string> bc)
  { for (int i=0; i<4; ++i) {bc[i][0] = tolower(bc[i][0]);} boundary_conditions = bc; }

  void print_boundary_conditions_to_screen()
  { for (int i=0; i<4; ++i) { cout << "bc["<<i<<"]: "<<boundary_conditions[i] << endl;} }

  /// @brief set the time step
  void set_timeStep( float dt )        { timeStep = dt; }

  /// @brief set the maximum time step
  void set_maxtimeStep (float max_dt)  { maxtimeStep = max_dt; }

  /// @brief set the ending time
  void set_endTime( float time )        { endTime = time; }

  /// @brief set the number of runs. Used for running multiple simulations from
  /// the same starting conditions
  void set_num_runs( int num )        { num_runs = num; }

  /// @brief sets the uplift mode
  void set_uplift_mode( int new_uplift_mode)    { uplift_mode = new_uplift_mode; }

  /// @brief overloaded function, set the array of uplift
  void set_uplift( Array2D <float> uplift )    { uplift_field = uplift; }

  /// @brief overloaded function, set the array of uplift, but using the uplift mode
  /// @details uplfit modes are:
  ///   (0) - block uplift
  ///  1   - tilt block
  ///  2   - gaussian
  ///  3   - quadratic
  /// 4   - periodic
  void set_uplift( int mode, float max_rate )
  { uplift_field = generate_uplift_field( mode, max_rate ); this->max_uplift = max_rate; }

  /// This adjusts the uplift mode
  void set_periodic_uplift(double uplift_amplitude_fraction)
          { uplift_amplitude = max_uplift*uplift_amplitude_fraction; uplift_mode = 4; }

  /// This sets the uplift amplitude as a fraction of the uplift
  void set_uplift_amplitude( double uplift_amplitude_fraction)
          { uplift_amplitude = max_uplift*uplift_amplitude_fraction; }

  /// @brief this sets the baseline uplift rate for the tilt block
  void set_baseline_uplift( float new_rate )    { baseline_uplift = new_rate; }

  /// @brief set the tolerance for determining steady state
  void set_steady_state_tolerance( float tol )    { steady_state_tolerance = tol; }

  /// @brief set the amplitude of random noise
  void set_noise( float noise_amp) { this->noise = noise_amp; }

  /// @brief sets fluvial erodibility
  void set_K( float K )          { this->K_fluv = K; }

  /// @brief sets the hillslope diffusivity
  void set_D( float D )          { this->K_soil = D; }

  /// @brief set the flexural rigidity
  void set_rigidity( float D )        { this->rigidity = D; }

  /// @brief sets the Area exponent in the SPIM
  void set_m( float m )          { this->m = m; }

  /// @brief sets the slope exponent in the SPIM
  void set_n( float n )          { this->n = n; }

  /// @brief sets the critical drainage area for channels
  void set_threshold_drainage( float area )    { this->threshold_drainage = area; }

  /// @brief Sets the critical slope
  void set_S_c( float Sc_new )        { S_c = Sc_new; }

  /// @brief Sets the periodicity in years
  void set_periodicity( float time )      { periodicity = time; }

  /// @brief Sets the 2nd periodicity in years
  void set_periodicity_2( float time )      { periodicity_2 = time; }

  /// @brief this snaps the periodicity to the timestep to ensure the max
  /// and min values of a varying parameter are reached
  void snap_periodicity( void );

  /// set the print interval
  void set_print_interval( int num_steps )    { print_interval = num_steps; float_print_interval = float(num_steps)*timeStep;  }

  /// set the float print interval
  //void set_float_print_interval( float float_dt_print )    { float_print_interval = float_dt_print; }

  /// set the float print interval
  void set_next_printing_time ( float next_float_dt_print )    { next_printing_time = next_float_dt_print; }

  /// @brief Update the raster data. WARNING this does not check the geometry of the raster
  /// @param Raster Another raster. It needs to be the same dimensions as the original raster
  /// @author SMM
  /// @date 16/12/2019
  void set_raster_data(LSDRaster& Raster);


  /// @brief makes a raster with a constant value
  /// @detail a brute force way to make K and U rasters when those are needed
  /// @param value the value that the raster pixels will take. 
  /// @author SMM
  /// @date 24/06/2021
  LSDRaster make_constant_raster(float value);


  /// @brief this sets the K mode
  /// @param mode The mode of calculating K
  /// K_mode == 1 sine wave
  /// K_mode == 2 square wave
  /// K_mode == 3 read from file
  /// K_mode: default is constant value
  /// @author JAJ
  /// @date 01/01/2014
  void set_K_mode( short mode )        { K_mode = mode; }

  /// @brief this sets the D mode
  /// @param mode The mode of calculating D
  /// D_mode == 1 sine wave
  /// D_mode == 2 square wave
  /// D_mode == 3 read from file
  /// D_mode: default is constant value
  /// @author JAJ
  /// @date 01/01/2014
  void set_D_mode( short mode )        { D_mode = mode; }

  /// @brief This sets the way the periodicity is calculated
  /// @param mode the mode of periodic forcing
  /// 1 (default) one periodicity used without
  /// 2 Two periodicities that switch at a given interval
  /// 3 Two periodicities used as a compound sin wave
  /// 4 Same as three, but weightings switch at a given interval (as in 2)
  /// @author JAJ
  /// @date 01/01/2014
  void set_period_mode( short mode )      { period_mode = mode; }

  /// @brief set the name of the model run
  void set_name( string name )        {this->name = name;}

  /// @brief set the name of the report
  void set_report_name( string name )      {report_name = name;}

  /// @brief set current frame, this is used for printing
  void set_current_frame(int new_frame)  { current_frame = new_frame; }

  /// @brief this just sets the initial steady state to true so that
  /// the periodic functions can be run from a starting DEM
  /// It also forces the time delay ans switch delay so periodic functions start straight away
  void force_initial_steady_state()   { initial_steady_state = true;
                                        time_delay = 0;
                                        switch_delay = 0;
                                      }

   // -------------------------------------------------------------------
  //@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  // Setters for turning on/off model components
  //@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  // -------------------------------------------------------------------
  /// @brief set the fluvial switch
  /// @param on_status a boolean, true if on, false if off
  /// @author JAJ
  /// @ date 01/01/2014
  void set_fluvial( bool on_status )      { fluvial = on_status; }

  /// @brief set the hillslop switch
  /// @param on_status a boolean, true if on, false if off
  /// @author JAJ
  /// @ date 01/01/2014
  void set_hillslope( bool on_status )      { hillslope = on_status; }

  /// @brief set the hillslop linear or nonlinear switch
  /// @param on_status a boolean, true if on, false if off
  /// @author SMM
  /// @ date 03/07/2014
  void set_nonlinear( bool on_status )      { nonlinear = on_status; }

  /// @brief set the isostacy switch
  /// @param on_status a boolean, true if on, false if off
  /// @author JAJ
  /// @ date 01/01/2014
  void set_isostasy( bool on_status )      { isostasy = on_status; }

  /// @brief set the flexure switch
  /// @param on_status a boolean, true if on, false if off
  /// @author JAJ
  /// @ date 01/01/2014
  void set_flexure( bool on_status )      { flexure = on_status; }

  /// @brief set the quiet switch
  /// @param on_status a boolean, true if on, false if off
  /// @author JAJ
  /// @ date 01/01/2014
  void set_quiet( bool on_status )      { quiet = on_status; }

  // PRINTING OF RASTERS
  /// @brief Sets the print elevation
  /// @param boolean; if true prints elevation
  /// @author SMM
  /// @date 09/08/2017
  void set_print_elevation( bool do_I_print_elevation )      { print_elevation = do_I_print_elevation; }

  /// @brief Sets the print hillshade
  /// @param boolean; if true prints hillshade
  /// @author SMM
  /// @date 09/08/2017
  void set_print_hillshade( bool do_I_print_hillshade )      { print_hillshade = do_I_print_hillshade; }

  /// @brief Sets the print erosion
  /// @param boolean; if true prints erosion
  /// @author SMM
  /// @date 09/08/2017
  void set_print_erosion( bool do_I_print_erosion )      { print_erosion = do_I_print_erosion; }

  /// @brief Sets the End time mode
  /// @param short 0,1,2,3,4
  /// @author BG
  /// @date 22/01/2019
  void set_endTime_mode(short edm){endTime_mode = edm;}



  // -------------------------------------------------------------------
  //@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  // Getter methods
  //@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
  // -------------------------------------------------------------------



  /// @brief gets the name of the model run from the data members
  /// @return the name of the model run
  /// @author JAJ
  /// @date 01/01/2014
  string get_name( void )          { return name; }

  /// @brief This function gets the fluvial erodability. It has a number of switches
  /// that determine how K is calcualted.
  /// K_mode == 1 sine wave
  /// K_mode == 2 square wave
  /// K_mode == 3 read from file
  /// K_mode: default is constant value
  /// @author JAJ
  /// @date 01/01/2014
  float get_K( void );

  /// @brief This function gets the soil transport coefficient erodability.
  /// It has a number of switches that determine how D is calculated.
  /// D_mode == 1 sine wave
  /// D_mode == 2 square wave
  /// D_mode == 3 read from file
  /// D_mode: default is constant value
  /// @author JAJ
  /// @date 01/01/2014
  float get_D( void );

  /// @brief Gets the area exponent
  /// @author SMM
  /// @date 10/08/2017
  float get_m( void )    { return m; }

  /// @brief Gets the slope exponent
  /// @author SMM
  /// @date 10/08/2017
  float get_n( void )    { return n; }


  /// @brief this gets the current_time
  float get_current_time( void )   { return current_time; }

  /// @brief this gets the endTime
  float get_endTime( void)         { return endTime; }

  /// @brief Gets the timestep
  float get_timeStep( void)        { return timeStep; }

  /// @brief Gets the print interval in years
  float get_float_print_interval( void)        { return float_print_interval; }

  /// @brief Gets the maximum timestep
  float get_maxtimeStep( void )    { return maxtimeStep; }

  /// @brief this gets the current frame for printing
  int get_current_frame( void)     { return current_frame;}

  /// @brief gets the uplift mode
  int get_uplift_mode(void)        { return uplift_mode; }


  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  // TOOLS FOR WRITING REPORTS AND PRINTING THINGS
  // @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief This method calculates some features of the landscape at set times. The
  /// frequency of the reports are set by the data member report_delay
  /// One of the things it does is calculates erosion rates, and stores this
  /// as a data member
  /// @author JAJ
  /// @date 01/01/2014
  void write_report( void );

  /// @brief this prints parameters to screen
  /// @author JAJ
  /// @date 01/01/2014
  void print_parameters( void );

  /// @brief this prints a file about what has happened over a cycle
  /// @author JAJ
  /// @date 01/01/2014
  void cycle_report( float, float, float);

  /// @brief this prints a final report (SMML not sure what is in the final report)
  /// @author JAJ
  /// @date 01/01/2014
  void final_report( void );

  /// @brief This function prints a series of rasters.
  /// The rasters printed depend on the switches
  /// print_elevation,  print_erosion, print_erosion_cycle, print_hillshade;
  /// and print_slope_area
  /// The filename inculdes the frame_num
  /// @param the frame of the rasters to be printed
  /// @author JAJ
  /// @date 01/01/2014
  void print_rasters( int frame_num );

  /// @brief This function prints a series of rasters.
  /// The rasters printed depend on the switches
  /// print_elevation,  print_erosion, print_erosion_cycle, print_hillshade;
  /// and print_slope_area
  /// The filename inculdes the frame_num
  /// It also prints a csv of the model info which can be ingested by pandas
  /// for visualisation
  /// @param the frame of the rasters to be printed
  /// @author FJC
  /// @date 22/08/17
  void print_rasters_and_csv( int frame );

  /// @brief This function prints the apparent cosmogenic rates from a collection
  /// of LSDParticle Columns
  /// @detail This function opens a file if none exists
  /// @param frame the frame to be printed
  /// @param CRNColumns the columns of cosmogenic particles
  /// @param  CRNParams and LSDCRNParamters object
  /// @author SMM
  /// @date 01/08/2014
  void print_average_erosion_and_apparent_erosion( int frame,
                                 vector<LSDParticleColumn>& CRNColumns,
                                 LSDCRNParameters& CRNParams);

  /// @brief This function prints the apparent cosmogenic rates from individual
  /// LSDParticle Columns
  /// @detail This function opens a file if none exists
  /// @param frame the frame to be printed
  /// @param CRNColumns the columns of cosmogenic particles
  /// @param  CRNParams and LSDCRNParamters object
  /// @author SMM
  /// @date 24/05/2015
  void print_column_erosion_and_apparent_erosion( int frame,
                                 vector<LSDParticleColumn>& CRNColumns,
                                 LSDCRNParameters& CRNParams);


  /// @brief This function closes some static outfiles used for printing
  /// @author SMM
  /// @date 01/08/2014
  void close_static_outfiles();

  /// @brief Print slope area data
  /// Probably fits better into LSDRaster, but requires LSDFlowInfo
  /// @param requires a filename
  /// @author JAJ
  /// @date 01/01/2014
  void slope_area_data( string name );

  /// @brief Print slope area data
  /// This is an overloaded function that calcualtes slope area
  /// data based on flags. There are two flag, one for the slope
  /// calculation and one for the area calculation
  /// @details The output is a file with four columns
  /// Elevation  slope area  predicted_slope
  /// The predicted slope is that based on the SS solution of stream power
  /// at the current uplift rate and the current K
  /// @param requires a filename
  /// @param a flag for calculation of the topographic slope
  /// 0 == polyfit using the data resolution as the smoothing diameter
  /// 1 == slope calculated with the D8 slopes, with dx = data resolution
  /// or DataResolution*sqrt(2) depending on flow direction.
  /// @param a flag for calculation of the area
  /// 0 == area using contributing pixels but smoothed to data resolution with polyfit
  /// 1 == area using contributing pixels only
  /// 2 ==
  /// @author SMM
  /// @date 18/06/2014
  void slope_area_data( string name, int slope_flag, int area_flag  );

  /// @brief Produce a template of a parameter file to be supplied to the model
  /// @param the name of the parameter file to be printed
  /// @author JAJ
  /// @date 01/01/2014
  void make_template_param_file(string filename);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // ~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~
  // MUDDPILE
  // nonlinear hillslope solver
  // Uses the boost library and mtl
  // ~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief this initiates some parameters for the assembler matrix
  /// it is required before running any further calculations for the
  /// nonlinear solver
  /// @details The function sets data members inv_dx_S_c_squared, dx_front_term,
  /// problem_dimension,  vec_k_value_i_j, vec_k_value_ip1_j,
  /// vec_k_value_im1_j, vec_k_value_i_jp1 and vec_k_value_i_jm1
  /// @author SMM
  /// @date 01/07/2014
  void MuddPILE_initiate_assembler_matrix( void );

  /// @brief this initiates some parameters for the assembler matrix
  /// it is required before running any further calculations for the
  /// nonlinear solver
  /// @details Does the work of calculating vec_k_value_i_j, vec_k_value_ip1_j,
  /// vec_k_value_im1_j, vec_k_value_i_jp1 and vec_k_value_i_jm1.
  /// These are the indices into the vectorized matrix of
  /// zeta values that are used in the assembly matrix
  /// the number of elements in the k vectors is N_rows*N_cols
  /// WARNING: This is only used for fixed NS boundaries and periodic EW boundaries
  /// @author SMM
  /// @date 01/07/2014
  void MuddPILE_calculate_k_values_for_assembly_matrix( void );

  /// @brief this assembles the sparse matrix that must then be solved
  /// to get the next iteration of the hillslope elevations
  /// @details Note this is not a buffered surface, boundaries are
  /// implemented at Row == 0 and Row = NRows-1 in the LSDModelRasters domain
  /// @param uplift_rate a float array giving the uplift rate
  /// @param fluvial_erosion_rate a float array giving the fluvial erosion rate
  /// @param mtl_Assembly_matrix this is a sparce matrix that is reset within
  /// this member function and passed to the solver
  /// @param mtl_b_vector the b vector in the linear system M z = b where
  /// M is the assembly matrix and z is the vector of surface elevations
  /// @author SMM
  /// @date 01/07/2014
  void MuddPILE_assemble_matrix(Array2D<float>& uplift_rate,
             Array2D<float>& fluvial_erosion_rate,
             mtl::compressed2D<float>& mtl_Assembly_matrix,
             mtl::dense_vector<float>& mtl_b_vector);

  /// @brief this function solves the assembled matrix for the nonlinear
  /// hillslope sediment flux law. The implementation calls
  /// MuddPILE_assemble_matrix.
  /// @details After solution the zeta_this_iter data member will be
  /// updated
  /// @param uplift_rate a float array of the same size as RasterData
  /// that contains the uplift rate
  /// @param fluvial_erosion_rate a float array of the same size as RasterData
  /// @author SMM
  /// @date 01/07/2014
  void MuddPILE_solve_assembler_matrix(Array2D<float>& uplift_rate,
             Array2D<float>& fluvial_erosion_rate);

  /// @brief This runs one timestep of the nonlinear sediment flux law
  /// It replaces the data in RasterData
  /// @param uplift_rate a float array of the same size as RasterData
  /// that contains the uplift rate
  /// @param fluvial_erosion_rate a float array of the same size as RasterData
  /// @param iteration_tolerance the maximum change in surface elevation
  /// between iterations
  /// @author SMM
  /// @date 01/07/2014
  void MuddPILE_nonlinear_creep_timestep(Array2D<float>& uplift_rate,
            Array2D<float>& fluvial_erosion_rate,
            float iteration_tolerance);

  /// @brief This version of the MuddPILE nonlinear solver does not include
  /// uplift and uses a default iteration tolerance of 1e-7. It is built
  /// to integrate with JAJ's 'run_components' module
  /// @author SMM
  /// @date 03/07/2014
  void MuddPILE_nl_soil_diffusion_nouplift();

  // some getters
  /// @return Number of rows as an integer.
  int get_NRows() const        { return NRows; }
  /// @return Number of columns as an integer.
  int get_NCols() const        { return NCols; }

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // Maps holding parameters and switches read from input paramter files
  //
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DAV - I am not sure whether these should be public or protected data members
  // LSDModelDriver needs read/write access to them, so my guess would be public?
  // DAV - Use getter/setter functions?

  // RM = LSDRasterModel specific
  // (CM = LSDCatchmentModel specific)

  /// This map holds all the possible model switches
  map<string,bool> RM_model_switches;

  /// This holds names of methods. For example, if the key is 'sed_transport_law', the string is
  /// the method which is used to calculate sediment transport (such as 'wilcock' or 'einstein')
  map<string,string> RM_method_map;

  /// This holds float parameters
  map<string,float> RM_float_parameters;

  /// This holds integer parameters
  map<string,int> RM_int_parameters;

  /// This holds names of supporting files, for example files that contain
  /// node of junction indices to be loaded.
  map<string,string> RM_support_file_names;



  protected:
  // Various parameters used in throughout the model run
  // These are set to the default values with default_parameters( void )
  // and then a parameter file is read to overwrite any of these explicitly set

  /// True if Supress output messages
  bool      quiet;

  ///  True if initialize_model has been run
  bool      initialized;

  /// True if initial steady state has been arrived at
  bool      steady_state;

  /// True for first steady state arrival, used for activating periodic forcing parameters
  bool      initial_steady_state;

  /// True for steady state using periodic forcing (not sure about this one)
  bool      cycle_steady_check;

  /// True if recording data
  bool      recording;

  /// True if writing the run report
  bool      reporting;

  /// Boundary conditions of model NESW
  vector <string>    boundary_conditions;

  /// Name of the model run
  string      name;

  /// Name of the report
  string       report_name;

  /// The current time
  float       current_time;

  /// The time at which steady state was reached. Used to sync period after steady state
  float      time_delay;

  /// Time in between each calculation
  float       timeStep;

  /// The maximum possible timestep. Used with adaptive timestepping
  float maxtimeStep;

  /// ending time
  float      endTime;

  /// This is a mode of end times. 0 == default, run until the endTime
  /// 1 == run until a specified time after steady state
  /// 2 == run a number of cycles
  /// 3 == run to steady state, then run some cycles.
  short      endTime_mode;

  /// the number of runs, used to do repeated cycles from the same steady state
  int      num_runs;

  /// The uplift field, can be used as absolute uplift or uplift rate, depending
  /// on the calling member function
  Array2D <float>   uplift_field;

  /// The mode of uplift. Options are:
  /// default == block uplift
  /// 1 == tilt block
  /// 2 == gaussian
  /// 3 == quadratic
  /// 4 == periodic
  int       uplift_mode;

  /// the maximum uplift rate
  float      max_uplift;

  /// the amplitude of uplift
  float uplift_amplitude;

  /// the baseline uplift rate (for different uplift modes)
  float baseline_uplift;

  /// an iteration tolerance for detemring if a model run is at steady state.
  float      steady_state_tolerance;
  float      steady_state_limit;

  /// Area exponent from the stream power law
  float     m;

  /// Slope exponent for the stream power law
  float     n;

  /// Fluvial erodability coefficient (units depend on m and n)
  float      K_fluv;

  /// Soil transport coefficient (usually in m^2/yr)
  float     K_soil;

  /// Drainage area above which soil will be flushed from the system in m^2
  float      threshold_drainage;

  /// Critical slope (for non-linear soil creep), dimensionless
  float      S_c;

  /// Flexural rigidity of plate (SMM: dimensions?)
  float      rigidity;

  /// Depth to topographic root, used in isostatic calculations
  Array2D <float>    root_depth;      // Depth of topographic root

  //  Measures of landscape response
  /// the erosion (distance)
  float      erosion;

  /// the erosion (distance) over the last timestep
  float      erosion_last_step;

  /// the erosion over cycles, a vector since it stores sucessive cycles
  vector<float>    erosion_cycle_record;

  /// Total erosion, calcualted as Offset from uplift from each cell
  float      total_erosion;

  /// minimum erosion distance, used in cyclic calculations
  float      min_erosion;

  /// minimum erosion distance, used in cyclic calculations
  float     max_erosion;

  /// maximum response over a single run (SMM: no idea if this is a length or what)
  float     response;

  /// response over all model runs (to be divided by num_runs)
  float      total_response;

  /// This sets the amplitude of random noise, used in model initialisation
  float      noise;

  /// SMM: not sure what this does
  float      report_delay;

  /// the elevations from the last timestep
  Array2D <float>    zeta_old;

  /// the elevations at steady state
  Array2D <float>    steady_state_data;

  /// a field calculated over an erosion cycle
  Array2D <float>    erosion_cycle_field;

  // Parameters for periodic forcing components
  /// Whether or not K is periodic
  /// K_mode == 1 sine wave
  /// K_mode == 2 square wave
  /// K_mode == 3 read from file
  /// K_mode: default is constant value
  short K_mode;

  /// Whether or not D is periodic
  /// D_mode == 1 sine wave
  /// D_mode == 2 square wave
  /// D_mode == 3 read from file
  /// D_mode: default is constant value
  short D_mode;

  /// Whether or not there is only one periodicity or two
  /// 1 (default) one periodicity used without
  /// 2 Two periodicities that switch at a given interval
  /// 3 Two periodicities used as a compound sin wave
  /// 4 Same as three, but weightings switch at a given interval (as in 2)
  short period_mode;

  /// Amplitude of K wave
  float K_amplitude;

  /// Amplitude of D wave
  float D_amplitude;

  /// Periodicty
  float periodicity;

  /// 2nd periodicity (if period_mode = 2)
  float periodicity_2;

  /// cycle that the model is on. Used to track changes in erosion between cycles
  int   cycle_number;

  /// Ratio for weight of periodicity in period_mode 3 or 4
  float p_weight;

  /// Time at which switch happens (time mode is same as endTime_mode)
  float switch_time;

  /// Similar to time delay (but for switching periodicities)
  float switch_delay;

  // Component switches
  /// True if fluvial erosion is on
  bool      fluvial;
  /// true if linear hillslope erosion is on
  bool      hillslope;
  /// True if nonlinear hillslope erosion is on
  bool      nonlinear;
  /// True if isostatic component is on
  bool      isostasy;
  /// True if // Whether flexural isostasy will be used
  bool      flexure;

  // Printing Utilities
  /// This is the current frame, used for keeping track of the output rasters
  int current_frame;

  /// interval over which output is written. Just based on number of timesteps
  int        print_interval;      // Interval at which output is written

  /// this is for printing at fixed times
  float float_print_interval;
  float next_printing_time;

  /// Switch for printing elevation, if true elevation is printed to a raster
  bool      print_elevation;

  /// Switch for printing erosion, if true elevation is printed to a raster
  bool      print_erosion;

  /// Switch for printing erosion over a cycle, if true erosion is printed to a raster
  bool      print_erosion_cycle;

  /// Switch for printing the hillshade, if true hillshade is printed to a raster
  bool      print_hillshade;

  /// Switch for printing S-A data, if true S-A data is printed
  bool      print_slope_area;

  /// This array keeps track of the elevation on the previous iteration
  Array2D<float> zeta_last_iter;

  /// This array keeps track of the elevation from the previous timestep
  Array2D<float> zeta_last_timestep;

  /// Array for the current iteratation (used for implicit nonlinear solvers)
  Array2D<float> zeta_this_iter;

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // Data members for the MuddPILE implicit solver
  //
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// precalculated 1/(dx^2*S_c^2)
  float inv_dx_S_c_squared;

  /// precalculated dt*D_nl/(dx^2)
  float dx_front_term;

  /// the problem dimension, used to determine the size of the solver matrix
  int problem_dimension;

  /// k value at i,j. This is an index into the vectorised elevation data.
  vector<int> vec_k_value_i_j;

  /// k value at i+1, j. This is an index into the vectorised elevation data.
  vector<int> vec_k_value_ip1_j;

  /// k value at i-1, j. This is an index into the vectorised elevation data.
  vector<int> vec_k_value_im1_j;

  /// k value at i, j+1. This is an index into the vectorised elevation data.
  vector<int> vec_k_value_i_jp1;

  /// k value at i, j-1. This is an index into the vectorised elevation data.
  vector<int> vec_k_value_i_jm1;



  private:
  void create();
  void create(string master_param);
  void create(string filename, string extension);
  void create(int ncols, int nrows, float xmin, float ymin,
        float cellsize, float ndv, Array2D<float> data, map<string,string>);
  void create(LSDRaster& An_LSDRaster);
  void default_parameters( void );


  /// @brief This function calculates the value of a sinusoidal periodic variable.
  /// To calcualte the variable, it uses the data member current_time
  /// and delay_switch to calcualte where in the cycle the parameter is.
  /// It uses a 'period mode' variable to determine the nature of the periodic
  /// forcing (which I'll have to look up later, SMM)
  /// @param the average value of the parameter
  /// @param the amplitude of variation in the parameter
  /// @return the value of the parameter at the current time
  /// @author JAJ
  /// @date 01/01/2014
  float periodic_parameter( float base_param, float amplitude );

  /// @brief This function calculates the value of a square wave periodic variable.
  /// To calcualte the variable, it uses the data member current_time
  /// and delay_switch to calcualte where in the cycle the parameter is.
  /// @param the average value of the parameter
  /// @param the amplitude of variation in the parameter
  /// @return the value of the parameter at the current time
  /// @author JAJ
  /// @date 01/01/2014
  float square_wave_parameter( float base_param, float amplitude );

  /// @brief load the K parameter from a stream
  /// The stream comes from a file called 'K_file' (with no extension)
  /// @return the K parameter at this timestesp
  /// @author JAJ
  /// @date 01/01/2014
  float stream_K_fluv(void);    // I don't like that these two are split up, but I was rushed and couldn't figure out how to pass an ifstream

  /// @brief load the D parameter from a stream
  /// The stream comes from a file called 'D_file' (with no extension)
  /// @return the D parameter at this timestesp
  /// @author JAJ
  /// @date 01/01/2014
  float stream_K_soil(void);

  /// @brief This is the main wrapper for running the landscape evolution model
  /// The nature of the model run is set by various switches (e.g., fluvial,
  /// hillslope, nonlinear hillslope, etc) and paramters that are contained
  /// within the object as data members.
  /// @author JAJ
  /// @date 01/01/2014
  void _run_components( void );
};

#endif
