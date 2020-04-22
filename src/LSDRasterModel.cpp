//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterModel.cpp
// cpp file for the LSDRasterModel object
// LSD stands for Land Surface Dynamics
// This object provides an environment for landscape evolution modelling, which can then
// be integrated with the topographic analysis tools to efficiently analyse model runs.
//
// The landscape evolution model uses implicit methods to provide stability with
// relatively long timesteps.  Fluvial erosion is solved following Braun and Willet (2013)
// using the fastscape algorithm, whilst hillslope sediment transport is modelled as a
// non-linear diffusive sediment flux, following the implicit scheme developed for
// MuDDPile.
//
// The aim is to have two complimentary models:
// i) a simple coupled hillslope-channel model in which large scale landscape dynamics
// can be modelled
// ii) a more complex treatment of hillslopes explicitly incorporating the role of
// vegetation in driving sediment production and transport, and that copes with the
// with the transition from soil mantled-bedrock hillslopes at high erosion rates.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
// Stuart Grieve, University of Edinburgh
// James Jenkinson, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1    24/07/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <string.h>
#include <queue>
#include <sys/stat.h>
#include <ctime>
#include <cstdlib>
#include <complex>
#include "TNT/tnt.h"
#include "TNT/jama_lu.h"
#include "TNT/jama_eig.h"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDRasterSpectral.hpp"
#include "LSDStatsTools.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDRasterModel.hpp"
#include "LSDSpatialCSVReader.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDParticleColumn.hpp"
#include "LSDJunctionNetwork.hpp"
using namespace std;
using namespace TNT;
using namespace JAMA;

#define PI 3.14159265358
#ifndef LSDRasterModel_CPP
#define LSDRasterModel_CPP

// The assignment operator
LSDRasterModel& LSDRasterModel::operator=(const LSDRasterModel& rhs)
 {
  if (&rhs != this)
   {
    create(rhs.get_NRows(),rhs.get_NCols(),rhs.get_XMinimum(),rhs.get_YMinimum(),
           rhs.get_DataResolution(),rhs.get_NoDataValue(),rhs.get_RasterData(), 
           rhs.get_GeoReferencingStrings());
   }
  return *this;
 }

// The destructor method, used to clean up temporary files and dynamic memory
// Currently only used to delete K/D files
LSDRasterModel::~LSDRasterModel( void )
{
  stringstream ss;
  if (K_mode == 3)
  {
    ss << ".K_file_" << name << ".aux";
    remove( ss.str().c_str() );
  }
  if (D_mode == 3)
  {
    ss << ".D_file_" << name << ".aux";
    remove( ss.str().c_str() );
  }

  close_static_outfiles();
}

// the create function.
// This sets up a model domain with a default size and model parameters
// Imposes UTM zone 1
void LSDRasterModel::create()
{
  NRows = 100;
  NCols = 100;
  DataResolution = 10;
  NoDataValue = -9999;
  XMinimum = 0;
  YMinimum = 0;
  RasterData = Array2D <float> (NRows, NCols, 0.0);

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);

  default_parameters();
}

// this creates a raster using an infile
void LSDRasterModel::create(string filename, string extension)
{
  read_raster(filename,extension);
}

// this creates a raster filled with no data values
void LSDRasterModel::create(int nrows, int ncols, float xmin, float ymin,
            float cellsize, float ndv, Array2D<float> data,
            map<string,string> GRS)
{
  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;
  GeoReferencingStrings =  GRS;

  RasterData = data.copy();

  if (RasterData.dim1() != NRows)
  {
    cout << "dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
  if (RasterData.dim2() != NCols)
  {
    cout << "dimension of data is not the same as stated in NCols!" << endl;
    exit(EXIT_FAILURE);
  }

  //int zone = 1;
  //string NorS = "N";
  //impose_georeferencing_UTM(zone, NorS);

}

// this creates a LSDRasterModel raster from another LSDRaster
void LSDRasterModel::create(LSDRaster& An_LSDRaster)
{
  NRows = An_LSDRaster.get_NRows();
  NCols = An_LSDRaster.get_NCols();
  XMinimum = An_LSDRaster.get_XMinimum();
  YMinimum = An_LSDRaster.get_YMinimum();
  DataResolution = An_LSDRaster.get_DataResolution();
  NoDataValue = An_LSDRaster.get_NoDataValue();
  GeoReferencingStrings =  An_LSDRaster.get_GeoReferencingStrings();
  RasterData = An_LSDRaster.get_RasterData();
}


LSDRasterModel::LSDRasterModel(int NRows, int NCols)
{
  this->NRows = NRows;
  this->NCols = NCols;
  this->DataResolution = 10;
  this->NoDataValue = -9999;
  XMinimum = 0;
  YMinimum = 0;
  RasterData = Array2D <float> (NRows, NCols, 0.0);

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);
}

// this creates an LSDRasterModel using a master parameter file
void LSDRasterModel::create(string master_param)
{
  NRows = 100;
  NCols = 100;
  DataResolution = 10;
  NoDataValue = -9999;
  XMinimum = 0;
  YMinimum = 0;

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);


  default_parameters();
  initialize_model(master_param);
}

// this sets default parameters for the model
void LSDRasterModel::default_parameters( void )
{
  initialized = false;
  name = "LSDRM";
  report_name = "LSDRM";
  reporting = true;
  vector <string> bc(4, "n");          // Initialise boundaries to No flow
  bc[0] = "b";
  bc[1] = "p";
  bc[2] = "b";
  bc[3] = "p";
  set_boundary_conditions( bc );          // Set these as default boundary conditions

  set_uplift( 0, 0.0005 );          // Block uplift, 0.005m.yr^{-1}
  set_baseline_uplift( 0.0005 );

  set_timeStep( 100 );            // 100 years
  set_maxtimeStep (1000);          // this limits the adaptive timestep
  set_endTime( 10000 );


  endTime_mode = 0;
  set_num_runs( 1 );
  set_K( 0.0002 );
  set_D( 0.02 );
  set_rigidity( 1E7 );
  set_m( 0.5 );
  set_n( 1 );
  set_threshold_drainage( -99 );          // Not used if negative
  set_S_c( 1.0 );              // 45 degrees, slope of 1

  set_print_interval( 10 );            // number of timesteps
  set_float_print_interval (5000);    // this is in years
  set_next_printing_time (0);

  set_steady_state_tolerance( 0.00001 );
  current_time = 0;
  noise = 0.1;

  K_mode = 0;
  D_mode = 0;
  periodicity   = 10000;
  periodicity_2 = 20000;
  period_mode = 1;
  switch_time = endTime/2;
  p_weight = 0.8;
  K_amplitude = 0.001;
  D_amplitude = 0.001;
  report_delay = 0;

  print_elevation = true;
  print_hillshade = false;
  print_erosion = false;
  print_erosion_cycle = false;
  print_slope_area = false;

  quiet =    false;
  fluvial =   true;
  hillslope =   true;
  nonlinear =  false;
  isostasy =   false;
  flexure =  false;

  steady_state_tolerance = 0.0001;
  steady_state_limit = -1;

  initialized = false;
  cycle_steady_check = false;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


// This returns the data in the raster model as a raster
LSDRaster LSDRasterModel::return_as_raster()
{
  LSDRaster NewRaster(NRows, NCols, XMinimum, YMinimum,
                      DataResolution, NoDataValue, RasterData,
                      GeoReferencingStrings);
  return NewRaster;
}

void LSDRasterModel::set_raster_data(LSDRaster& Raster)
{
  Array2D<float> temp_data = Raster.get_RasterData();
  RasterData = temp_data.copy(); 
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This adds a path to the run name and the report name
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::add_path_to_names( string pathname)
{
  string lchar = pathname.substr(pathname.length()-2,1);
  string slash = "/";
  cout << "lchar is " << lchar << " and slash is " << slash << endl;

  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;
    pathname = pathname+slash;
  }

  name = pathname+name;
  report_name = pathname+name;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// INITIALISATION MODULE
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This module initialises the model runs, calling the required function from
// the initial topography and loads the parameters from the parameter file.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  Initialise_model
//----------------------------------------------------------------------------
void LSDRasterModel::initialize_model( string& parameter_file, string& run_name,
  float& dt, float& EndTime, float& PrintInterval,
  float& k_w, float& b, float& m, float& n, float& K, float& ErosionThreshold,
  float& K_nl, float& S_c, float& UpliftRate, float& PrecipitationRate,
  float& NorthBoundaryElevation, float& SouthBoundaryElevation,
  Array2D<float>& PrecipitationFlux, Array2D<float>& SlopesBetweenRows,
  Array2D<float>& SlopesBetweenColumns, Array2D<float>& ErosionRate)
{
    // load the parameters
  // each parameter in the param file is preceded by its name
  // these MUST be in the correct order
  // the names MUST NOT have spaces
  string temp;
  ifstream param_in;
  param_in.open(parameter_file.c_str());

  param_in >> temp >> run_name;
  cout << "run name is: " << run_name;
  param_in >> temp >> dt;
  param_in >> temp >> EndTime;
  param_in >> temp >> PrintInterval;
  cout << "dt: " << dt << " end_time: " << EndTime << " print_interval: " << PrintInterval << endl;
  param_in >> temp >> k_w;
  param_in >> temp >> b;
  param_in >> temp >> m;
  param_in >> temp >> n;
  param_in >> temp >> K;
  param_in >> temp >> ErosionThreshold;
  cout << "k_w: " << k_w << " b: " << b << " m: " << m << " n: " << n << " K: " << K << " eros_thresh: " << ErosionThreshold << endl;
  param_in >> temp >> K_nl;
  param_in >> temp >> S_c;
  cout << "D_nl: " << K_nl << " S_c: " << S_c << endl;
  param_in >> temp >> UpliftRate;
  param_in >> temp >> PrecipitationRate;
  cout << "uplift_rate: " << UpliftRate << " precip_rate: " << PrecipitationRate << endl;
  param_in >> temp >> NorthBoundaryElevation;
  param_in >> temp >> SouthBoundaryElevation;
  cout << "N bdry elev: " << NorthBoundaryElevation << " S bdry elev: " << SouthBoundaryElevation << endl;

//   string surface_file;
//   param_in >> surface_file;
//
//   string file_extension;
//   param_in >> file_extension;
//
//   cout << "Surface_file is: " << surface_file << endl;

  param_in.close();

  float dx = get_DataResolution();
  float dy = get_DataResolution();
  cout << " NRows: " << NRows << " NCols: " << NCols << " dx: " << dx << " dy: " << dy
       << " xllcorn: " << XMinimum << " yllcorn: " << YMinimum << endl;

  // now set up some arrays
  // first the precipitation array
  PrecipitationFlux = precip_array_from_precip_rate(PrecipitationRate);

  // set up slope arrays of the correct size
  Array2D<float> slopes_between_rows_temp(NRows+1,NCols,0.0);
  Array2D<float> slopes_between_columns_temp(NRows,NCols+1,0.0);
  SlopesBetweenRows = slopes_between_rows_temp.copy();
  SlopesBetweenColumns = slopes_between_columns_temp.copy();

  // set up erosion rate array of the correct size
  Array2D<float> temp_erate(NRows,NCols,0.0);
  ErosionRate = temp_erate.copy();
}

// -----------------------------------------------------------------------
// Alternative mode of initialising from parameter file
//
// Loads parameters using the void parse_line method found in LSDStatsTools
// Parameters are loaded into intrinsic class attributes
// Part of a move to using attributes vs passing through functions
//  ----------------------------------------------------------------------
void LSDRasterModel::initialize_model(string param_file)
{
  bool loaded_from_file = false;
  initialized = true;
  ifstream infile;
  infile.open(param_file.c_str());
  string parameter, value, lower;

  while (infile.good())
  {
    parse_line(infile, parameter, value);
    lower = parameter;
    if (parameter == "NULL")
      continue;
    for (unsigned int i=0; i<parameter.length(); ++i)
      lower[i] = tolower(parameter[i]);

    if   (lower == "run name")    name     = value;
    else if (lower == "time step")    timeStep   = atof(value.c_str());
    else if (lower == "end time")    endTime   = atof(value.c_str());
    else if (lower == "num runs")    num_runs  = atoi(value.c_str());
    else if (lower == "end time mode")  endTime_mode  = atoi(value.c_str());
    else if (lower == "max uplift")    max_uplift = atof(value.c_str());
    else if (lower == "baseline uplift")  baseline_uplift = atof(value.c_str());
    else if (lower == "uplift mode")  uplift_mode   = atoi(value.c_str());
    else if (lower == "tolerance")    steady_state_tolerance = atof(value.c_str());
    else if (lower == "steady limit")  steady_state_limit = atof(value.c_str());
    else if (lower == "boundary code")  for (int i=0; i<4; ++i) boundary_conditions[i] = value[i];
    else if (lower == "m")      m     = atof(value.c_str());
    else if (lower == "n")      n     = atof(value.c_str());
    else if (lower == "k")      K_fluv     = atof(value.c_str());
    else if (lower == "threshold drainage") threshold_drainage = atof(value.c_str());
    else if (lower == "d")      K_soil     = atof(value.c_str());
    else if (lower == "s_c")    S_c     = atof(value.c_str());
    else if (lower == "rigidity")    rigidity  = atof(value.c_str());
    else if (lower == "nrows"){    if (loaded_from_file == false)   NRows     = atoi(value.c_str());}
    else if (lower == "ncols"){    if (loaded_from_file == false)   NCols     = atoi(value.c_str());}
    else if (lower == "resolution"){  if (loaded_from_file == false)   DataResolution   = atof(value.c_str()); }
    else if (lower == "print interval")  print_interval  = atoi(value.c_str());
    else if (lower == "k mode")    K_mode    = atoi(value.c_str());
    else if (lower == "d mode")    D_mode     = atoi(value.c_str());
    else if (lower == "periodicity")  periodicity   = atof(value.c_str());
    else if (lower == "periodicity 2")  periodicity_2   = atof(value.c_str());
    else if (lower == "P ratio")    {p_weight  = atof(value.c_str()); if (p_weight > 1) p_weight = 1;}
    else if (lower == "period mode")  period_mode  = atoi(value.c_str());
    else if (lower == "switch time")  switch_time  = atof(value.c_str());
    else if (lower == "k amplitude")  K_amplitude  = atof(value.c_str()) * K_fluv;
    else if (lower == "d amplitude")  D_amplitude   = atof(value.c_str()) * K_soil;
    else if (lower == "noise")    noise    = atof(value.c_str());
    else if (lower == "report delay")  report_delay   = atof(value.c_str());
    else if (lower == "fluvial")    fluvial   = (value == "on") ? true : false;
    else if (lower == "hillslope")    hillslope   = (value == "on") ? true : false;
    else if (lower == "non-linear")    nonlinear   = (value == "on") ? true : false;
    else if (lower == "isostasy")    isostasy   = (value == "on") ? true : false;
    else if (lower == "flexure")    flexure   = (value == "on") ? true : false;
    else if (lower == "quiet")    quiet    = (value == "on") ? true : false;
    else if (lower == "reporting")    reporting  = (value == "on") ? true : false;
    else if (lower == "print elevation")  print_elevation = (value == "on") ? true : false;
    else if (lower == "print hillshade")  print_hillshade = (value == "on") ? true : false;
    else if (lower == "print erosion")
    {
      print_erosion   = (value == "on") ? true : false;
      cout << "I will print the erosion!" << endl;
    }
    else if (lower == "print erosion cycle") print_erosion_cycle = (value == "on") ? true : false;
    else if (lower == "print slope-area")  print_slope_area= (value == "on") ? true : false;

    else if  (lower == "load file")
    {
      ifstream file(value.c_str());
      if (file)
      {
        file.close();
        read_raster(value.substr(0, value.find(".")), value.substr(value.find(".")+1));
        loaded_from_file = true;
      }
      else
        cerr << "Warning, file '" << value << "' not found" << endl;
    }

    else  cout << "Line " << __LINE__ << ": No parameter '" << parameter << "' expected.\n\t> Check spelling." << endl;
  }
  //if (hillslope)
  //  steady_state_tolerance *= pow(10, 2.8);
  if (name != "")
    report_name = name;
  else
    report_name = param_file;
  if (loaded_from_file == false)
  {
    RasterData = Array2D<float>(NRows, NCols, 0.0);
    // Generate random noise
    random_surface_noise(0, noise);
    // Fill the topography
    LSDRaster *temp;
    temp = new LSDRaster(*this);
    float thresh_slope = 0.00001;
    *temp = fill(thresh_slope);
    RasterData = temp->get_RasterData();
    delete temp;
  }
  root_depth = Array2D<float>(NRows, NCols, 0.0);
  current_time = 0;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function appends the run name, so that if you want you can add some
// details to the filename
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::append_run_name(string append_name)
{
  name = name+append_name;
  cout << "The model run name is: " << name << endl;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function adds a random float to each pixel in the raster
// Written sometime 2014 JAJ
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::random_surface_noise( float min, float max )
{

  cout << "Seeding the surface with random asperities of the amplitude " << max-min << endl;

  // Seed random numbers
  short dimension;
  int size;
  bool periodic;
  interpret_boundary(dimension, periodic, size);
  int start_i, end_i;
  int start_j, end_j;

  if (dimension == 0)
  {
    start_i = 1; end_i = NRows-2;
    start_j = 0; end_j = NCols-1;
  }
  else
  {
    start_i = 0; end_i = NRows-1;
    start_j = 1; end_j = NCols-2;
  }

  srand( static_cast <unsigned> (time(0)) );

  // Add random float to each pixel
  for (int i=start_i; i<=end_i; ++i)
  {
    for (int j=start_j; j<=end_j; ++j)
    {
      if (is_base_level(i, j))
        continue;
      RasterData[i][j] += static_cast <float> ( rand()) / static_cast <float> (RAND_MAX/(max-min)) + min;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Similar to above, but uses the noise parameter stored in the object
// SMM 17/6/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::random_surface_noise()
{
  // check on the noise data member. It needs to be bigger than 10-6)
  // if not set 1 mm as default
  if (noise < 0.000001)
  {
    noise = 0.001;
  }

  // we set the min and max between zero and noise. We don't go between
  // -noise/2 and noise/2 just because we don't want negative elevations near
  // a base level node.
  float min = 0;
  float max = noise;

  cout << "Seeding the surface with random asperities of the amplitude " << noise << endl;

  // Seed random numbers
  short dimension;
  int size;
  bool periodic;
  interpret_boundary(dimension, periodic, size);
  int start_i, end_i;
  int start_j, end_j;

  if (dimension == 0)
  {
    start_i = 1; end_i = NRows-2;
    start_j = 0; end_j = NCols-1;
  }
  else
  {
    start_i = 0; end_i = NRows-1;
    start_j = 1; end_j = NCols-2;
  }

  srand( static_cast <unsigned> (time(0)) );

  // Add random float to each pixel
  for (int i=start_i; i<=end_i; ++i)
  {
    for (int j=start_j; j<=end_j; ++j)
    {
      if (is_base_level(i, j))
        continue;
      RasterData[i][j] += static_cast <float> ( rand()) / static_cast <float> (RAND_MAX/(max-min)) + min;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// initializes a parabolic surface with elevations on the north and south edges at zero and
// elevation in the middle of 'peak elevation'
// the parabola also has random noise on it, with amplitude stored in
// the data member 'noise'
// The noise only adds to the elevations since we don't want elevations less
// than zero.
// Default noise is 1mm
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::initialise_parabolic_surface(float peak_elev, float edge_offset)
{
  // check on the noise data member. It needs to be bigger than 10-6)
  // if not set 1 mm as default
  if (noise < 0.000001)
  {
    noise = 0.001;
  }

  // set up the length coordinate
  float local_x;
  float L = DataResolution*(NRows-1);
  float row_elev;
  float perturb;

  // the seed for the random perturbation
  long seed = time(NULL);

  // loop through getting the parabolic elevation at each row, and then
  // writing across the entire domain
  for(int row = 0; row < NRows; row++)
  {

    local_x = row*DataResolution;
    row_elev = - 4.0*(local_x*local_x-local_x*L)*peak_elev / (L*L);

    for (int col = 0; col < NCols; col++)
    {

      // at N and S boundaries, the elevation is set to 0
      if( row == 0 || row == NRows-1)
      {
         RasterData[row][col] = 0;
      }
      else      // elsewhere initiate with a parabola
      {
        perturb = (ran3(&seed))*noise;
        RasterData[row][col] = row_elev + perturb + edge_offset;
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// superimposes a parabolic surface with elevations on the north and south edges at zero and
// elevation in the middle of 'peak elevation'
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::superimpose_parabolic_surface(float peak_elev)
{

  // set up the length coordinate
  float local_x;
  float L = DataResolution*(NRows-1);
  float row_elev;

  // loop through getting the parabolic elevation at each row, and then
  // writing across the entire domain
  for(int row = 0; row < NRows; row++)
  {

    local_x = row*DataResolution;
    row_elev = - 4.0*(local_x*local_x-local_x*L)*peak_elev / (L*L);

    for (int col = 0; col < NCols; col++)
    {

      // at N and S boundaries, tthere is no perturbation
      if( row != 0 && row != NRows-1)
      {
        RasterData[row][col] = RasterData[row][col]+row_elev;
      }
    }
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Initialies a fractal-based terrain surface with given fractal dimension, D
// D should be between 2 and 3 for 'reasonable' landscapes
// This appraoch is known as the Midpoint method
// More algorithms can be found in Saupe (1987): Algorithms for random fractals
// (an interesting bedtime read!)
//
// DAV 15/10/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//void LSDRasterModel::intialise_MFD_fractal_surface(float fractal_D)
//{
  //// set up the length coordinate
  //float local_x;
  //float L = DataResolution*(NRows-1);
  //float row_elev;
  //float perturb;

  //// Arguments for the Fractalombulator
  //float maxlevel;    // maximum number of recursions to use; N = 2^maxlevel
  //float sigma;      // intial standard deviation
  //float H;        // H is the parameter determining the fractal dimension; D = 3 - H
  //bool addition;      // boolean paramter that turns random additions on or off
  //float seed;        // seed for the random number gen

  //// Variables
  //int i, N, stage;   // i is just a wee counter deely. N is number of recursions, stage is the
  //float delta;    // holds the standard deviation
  //int x, y, y0, D, d;   // integers for the indexing of arrays

  ////Gaussian functions defined in LSDStatsTools

  //// Lambda function (new in C++11) to return the gauss random number
  //// the bit inside the []  are variables we want to 'capture' from the outer function.
  //// the bit inside the { } is the lambda function itself

  //double f3 = [delta, x0, x1, x2] { return ((x0 + x1 + x2)/3 + delta * Gauss_rand(100, 0.0, 1.0) }
  //double f4 = [delta, x0, x1, x2, x3] { return ((x0 + x1 + x2 + x3)/4 + delta * Gauss_rand(100, 0.0, 1.0) }


  ////Fractalombulator

  //N = pow(2, maxlevel);
  //delta = sigma;

  //Array2D<double> X(N+1,N+1);

  //// set the four corners to random numbers
  //X[0][0] = delta * Gauss_rand(100, 0.0, 1.0);    // this is a bad way to do it, set them as named global? variables
  //X[0][N] = delta * Gauss_rand(100, 0.0, 1.0);
  //X[N][0] = delta * Gauss_rand(100, 0.0, 1.0);
  //X[N][N] = delta * Gauss_rand(100, 0.0, 1.0);

  //D = N;
  //d = N/2;

  //for (stage = 1, stage <= maxlevel, stage++)
  //{
    //delta = delta * pow(0.5, 0.5*H)    // switch to the 45 degree grid type (See Saupe '87)

    //for (x=d; x <= (N-d); x++)
    //{
      //for (y=d; y <= (N-d); y++)
      //{
        //X[x][y] = f4(delta...    // TO COMPLETE
      //}
    //}
  //}


//}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This creates a fractal surface DEM using the Fourier filtering method
// (Saupe 1987)
//
// DAV 17/10/2014
//
// SMM 10/08/2017 NOT WORKING
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::intialise_fourier_fractal_surface(float fractal_D)
//
// N is size of the array along 1 dimension (NRows)
//
// A[][] is a 2D array of complex variables, size N^2
{

  int N = NRows;
  //int N = RasterData.get_NRows;
  float H = fractal_D;
  Array2D< std::complex<float> > A;
  Array2D<float> X;
  std::complex<float> this_complex;


  //int i0, j0;

  float rad, phase;  //Polar coordinates of the Fourier coefficient

  for (int i = 0; i<=(N/2); i++)
  {
    for (int j = 0; j<=(N/2); j++)
    {
      phase = 2 * PI * rand();

      if (i != 0 || j != 0)
      {
        rad = pow(i*i + j*j, (-(H+1)/2)) * Gauss_rand(100, 0.0, 1.0);
      }
      else
      {
        rad = 0;
      }

      this_complex.real(rad * cos(phase));
      this_complex.imag(rad * sin(phase));
      A[i][j] = this_complex;

      //A[i][j] = {rad * cos(phase), rad * sin(phase)};     // left of the comma is real part, right of the comma is imaginary part
      // This { } notation may only work in C++11 (i.e. the most recent standard as of 2014)



      /* This stuff seems to have no effect
      if (i==0)
      {
        i0 = 0;
      }
      else
      {
        i0 = N-i;
      }

      if (j==0)
      {
        j0 = 0;
      }
      else
      {
        j0 = N-j;
      }
      */

      this_complex.real(rad * cos(phase));
      this_complex.imag(-rad * sin(phase));
      A[i][j] = this_complex;

      // A[i0][j0] = {rad * cos(phase), -rad * sin(phase)};

    }
  }

  // Now for the *imaginary* numbers part
  // We are setting the imaginary parts of the midpoint-edges of the array to zero
  A[N/2][0].imag(0.0);
  A[0][N/2].imag(0.0);
  A[N/2][N/2].imag(0.0);

  // Now generate the fractal surface
  for (int i=1; i<=(N/2)-1; i++)
  {
    for (int j=1; j<=(N/2)-1; j++)
    {
      phase = 2 * PI * rand();
      rad = pow(i*i + j*j, -(H+1)/2) * Gauss_rand(100, 0.0, 1.0);

      this_complex.real(rad * cos(phase));
      this_complex.imag(rad * sin(phase));
      A[i][N-j] = this_complex;

      this_complex.real(rad * cos(phase));
      this_complex.imag(-rad * sin(phase));
      A[N-i][j] = this_complex;


      //A[i][N-j] = {rad * cos(phase), rad * sin(phase)};    // left of the comma is real part, right of the comma is imaginary part
      //A[N-i][j] = {rad * cos(phase), -rad * sin(phase)};

    }
  }

  // Here is the place to do the fast fourier transform: DAV
  // Note: There looks to be a very similar function in LSDRasterSpectral!
  // InvFFT2D(A,X,N)

  dfftw2D_inv_complex(A, X, 1);

  RasterData = X.copy();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This creates a fractal surface DEM using the Fourier filtering method
// (Saupe 1987)
//
// Like above, but uses the version in LSDRasterSpectral
// ONLY WORKS ON SMALL DEMS WITH pow(2) dimensions
//
// SMM 10/08/2017

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::intialise_fourier_fractal_surface_v2(float beta, float desired_relief)
{
  Array2D<float> zeta=RasterData.copy();

  // Create a raster spectral
  LSDRasterSpectral Fourier(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);

  Fourier.generate_fractal_surface_spectral_method(beta,desired_relief);

  zeta = Fourier.get_RasterData();

  RasterData = zeta.copy();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This creates a fractal surface DEM using the diamond square algorithm.
// Believe it or not I lifted this algorithm from Notch, the creator of Minecraft,
// who posted it online and then had it modified by Charles Randall
// https://www.bluh.org/code-the-diamond-square-algorithm/
//
// SMM 10/08/2017

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::intialise_diamond_square_fractal_surface(int feature_order, float desired_relief)
{
  Array2D<float> zeta=RasterData.copy();

  // temporary raster for performing diamond square
  LSDRaster temp_raster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);

  // get the dimaond square raster
  // IMPORTANT: this will be bigger than the original raster
  LSDRaster DSRaster = temp_raster.DiamondSquare(feature_order, desired_relief);

  // resample the raster to get a surface the correct size
  // it won't wrap but the running to steady will take care of that.
  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      zeta[row][col] = DSRaster.get_data_element(row,col);
    }
  }

  RasterData = zeta.copy();

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This tapers the north and south boundaries to 0 elevation and raises the
// entire DEM above sea level
//
// SMM 10/08/2017

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::initialise_taper_edges_and_raise_raster(int rows_to_taper)
{
  Array2D<float> zeta=RasterData.copy();

  // first we need to loop through all the data and raise above the sea level,
  // or lower to sea level accordingly
  float MinElev = 9999999;
  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if (zeta[row][col] < MinElev)
      {
        MinElev = zeta[row][col];
      }
    }
  }
  cout << "Tapering. Found the mininum elevation, it is: " << MinElev << endl;
  if (MinElev == -9999)
  {
    cout << "Your min elev is no data value. Setting it to 0. This means I won't raise or lower your raster" << endl;
    MinElev = 0;
  }


  // now adjust elevations so the lowst points are at zero elevation
  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      zeta[row][col] = zeta[row][col]-MinElev;
    }
  }

  // now we loop through the edge nodes, muliplying each by a fraction so they taper to zero elevation
  float this_frac;
  for (int taper_row = 0; taper_row < rows_to_taper; taper_row++)
  {
    this_frac = float(taper_row)/float(rows_to_taper);

    for(int col = 0; col<NCols; col++)
    {
      zeta[taper_row][col] = this_frac*zeta[taper_row][col];
      zeta[NRows-1-taper_row][col] = this_frac*zeta[NRows-1-taper_row][col];
    }
  }

  RasterData = zeta.copy();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This raises and then fills the DEM
//
// SMM 25/08/2017

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::raise_and_fill_raster()
{
  Array2D<float> zeta=RasterData.copy();

  // first we need to loop through all the data and raise above the sea level,
  // or lower to sea level accordingly
  float MinElev = 9999999;
  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if (zeta[row][col] < MinElev)
      {
        MinElev = zeta[row][col];
      }
    }
  }
  cout << "Raising raster. Found the mininum elevation, it is: " << MinElev << endl;


  // now adjust elevations so the lowst points are at zero elevation
  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      zeta[row][col] = zeta[row][col]-MinElev;
    }
  }

  RasterData = zeta.copy();

  cout << "Now I am filling the data" << endl;
  LSDRaster *temp;
  temp = new LSDRaster(*this);
  float thresh_slope = 0.00001;
  *temp = fill(thresh_slope);
  RasterData = temp->get_RasterData();
  delete temp;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This raises and then fills the DEM
// overloaded to take the min slope as an argument - FJC July 2018
//
// SMM 25/08/2017

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::raise_and_fill_raster(float min_slope_for_fill)
{
  Array2D<float> zeta=RasterData.copy();

  // first we need to loop through all the data and raise above the sea level,
  // or lower to sea level accordingly
  float MinElev = 9999999;
  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if (zeta[row][col] < MinElev)
      {
        MinElev = zeta[row][col];
      }
    }
  }
  cout << "Raising raster. Found the mininum elevation, it is: " << MinElev << endl;


  // now adjust elevations so the lowst points are at zero elevation
  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      zeta[row][col] = zeta[row][col]-MinElev;
    }
  }

  RasterData = zeta.copy();

  cout << "Now I am filling the data" << endl;
  LSDRaster *temp;
  temp = new LSDRaster(*this);
  *temp = fill(min_slope_for_fill);
  RasterData = temp->get_RasterData();
  delete temp;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This raises part of the raster by a specified amount, in metres. Simulates
// fault cutting horizontally across the raster.
// throw amount is instantaneously applied.
// int throw_type: 0 = top third of raster
//                 1 = top half of raster
//                 2 = top two thirds of raster
//
// FJC 06/07/18

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::normal_fault_part_of_raster(int throw_amt, int throw_type)
{
  Array2D<float> zeta=RasterData.copy();

  if (throw_type == 0)
  {
    for (int row = 0; row < int(NRows/3); row++)
    {
      for (int col = 0; col < NCols; col++)
      {
        zeta[row][col] = zeta[row][col] + throw_amt;
      }
    }
  }
  else if (throw_type == 1)
  {
    for (int row = 0; row < int(NRows/2); row++)
    {
      for (int col = 0; col < NCols; col++)
      {
        zeta[row][col] = zeta[row][col] + throw_amt;
      }
    }
  }
  else if (throw_type == 2)
  {
    for (int row = 0; row < int(NRows - NRows/3); row++)
    {
      for (int col = 0; col < NCols; col++)
      {
        zeta[row][col] = zeta[row][col] + throw_amt;
      }
    }
  }

  RasterData = zeta.copy();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This simulates a base level fall by raising the entire raster except the
// first and last rows by a certain number of metres
//
// FJC 18/07/18

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::base_level_fall(int uplift_amt)
{
  Array2D<float> zeta=RasterData.copy();

  for (int row = 0; row < NRows - 1; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      bool base_level = is_base_level(row, col);
      if (base_level == false)
      {
        zeta[row][col] = zeta[row][col] + uplift_amt;
      }
    }
  }

  RasterData = zeta.copy();
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::AdjustElevation(float elevation_change)
{
  for(int row = 0; row< NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        RasterData[row][col] = RasterData[row][col]+elevation_change;
      }
    }
  } 
} 


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this creates a hillslope at steady state for the nonlinear sediment flux law
// Solution from Roering et al., (EPSL, 2007)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::initialise_nonlinear_SS(float U)
{

  // check to make sure there are values for K_soil and S_c that are reasonable
  if (K_soil > 1 || K_soil < 1e-6)
  {
    cout << "Warning, LSDRasterModel::initialise_nonlinear_SS, D does not appear to be set" << endl;
    cout << "Defaulting to 0.0001 m^2/yr" << endl;
    K_soil = 0.0002;
  }
  if (S_c > 1.5 || S_c  < 0.3)
  {
    cout << "Warning, LSDRasterModel::initialise_nonlinear_SS, S_c does not appear to be set" << endl;
    cout << "Defaulting to 1 m^2/yr" << endl;
    S_c = 1.0;
  }

  float loc_y;
   float rho_ratio = 1;                    // need to double check if density conversion
                                           // is factored into the nonlinear solver (SMM, 2/7/2014)
   float hillslope_length = (NRows-1)*DataResolution;
   float divide_loc = hillslope_length*0.5;

  cout << "LSDRasterModel::initialise_nonlinear_SS, D_nl is: "
       << K_soil << " and S_c is: " << S_c << endl;

   cout << "Hillslope length is: " << hillslope_length
       << " and the divide location is: " << divide_loc << endl;

  // some precalculations to save time
  float beta = rho_ratio*U;
   float leading_term = -S_c*S_c/(2*beta);
   float inside_sqrt;
   float sqrt_term;
   float log_term;

  vector<float> elevs;
  for (int row = 0; row<NRows; row++)
  {
    // get the position along the slope
    loc_y = DataResolution*float(row);
    inside_sqrt = K_soil*K_soil+ (4*(loc_y-divide_loc)*(loc_y-divide_loc)*beta*beta/(S_c*S_c));
    sqrt_term = sqrt(inside_sqrt);
    log_term = log(((sqrt_term+K_soil)*S_c)/(2*beta) );
    elevs.push_back(leading_term*(sqrt_term-K_soil*log_term));

    cout << "Location is " << loc_y << " and elevation is " << elevs[row]-elevs[0] << endl;
  }


  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {

      RasterData[row][col] = elevs[row]-elevs[0];
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This resizes and resets the model
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::resize_and_reset( int new_rows, int new_cols )
{
  // set up some empty arrays
  Array2D<float> empty_array;
  Array2D<float> empty_array_sized(new_rows,new_cols,0.0);

  // set most of the arrays as data members to empty arrays
  uplift_field = empty_array.copy();
  root_depth =  empty_array.copy();
  zeta_old = empty_array.copy();
  steady_state_data = empty_array.copy();
  erosion_cycle_field = empty_array.copy();
  zeta_last_iter = empty_array.copy();
  zeta_last_timestep = empty_array.copy();
  zeta_this_iter = empty_array.copy();

  // reset the size of the RasterData
  RasterData = empty_array_sized.copy();

  // reset the rows and columns
  NRows = new_rows;
  NCols = new_cols;

  // now generate a random surface
  random_surface_noise();

  // reset flags and total erosion rates
  total_erosion = 0;
  total_response = 0;
  steady_state = false;
  initial_steady_state = false;

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This resizes and resets the model
// This overloaded version also resets the data resolution
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::resize_and_reset( int new_rows, int new_cols, float new_resolution )
{
  // set up some empty arrays
  Array2D<float> empty_array;
  Array2D<float> empty_array_sized(new_rows,new_cols,0.0);

  // set most of the arrays as data members to empty arrays
  uplift_field = empty_array.copy();
  root_depth =  empty_array.copy();
  zeta_old = empty_array.copy();
  steady_state_data = empty_array.copy();
  erosion_cycle_field = empty_array.copy();
  zeta_last_iter = empty_array.copy();
  zeta_last_timestep = empty_array.copy();
  zeta_this_iter = empty_array.copy();

  // reset the size of the RasterData
  RasterData = empty_array_sized.copy();

  // reset the rows and columns
  NRows = new_rows;
  NCols = new_cols;

  DataResolution = new_resolution;

  // now generate a random surface
  random_surface_noise();

  // reset flags and total erosion rates
  total_erosion = 0;
  total_response = 0;
  steady_state = false;
  initial_steady_state = false;

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
// FUNCTIONS FOR CHECKING STATE VARIABLES
//@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function checks to see if the model has arrived at steady state
// Written by JAJ sometime 2014. Commented by SMM 27/06/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::check_steady_state( void )
{
  // start by assuming steady_state has been acheived and then
  // reset flag to false if not.
  steady_state = true;

  //cout << "Line 591, checking steady state, at this point ISS is" << initial_steady_state << endl;

  //  SMM: this check if steady state has been arrived across a cycle.
  // It check a data member erosion_cycle_record, which is calculated elsewhere
  if (cycle_steady_check)
  {
    for (int i=0; i<4; ++i)
    {
      if (erosion_cycle_record[i] == -99 || abs(erosion_cycle_record[i]
             - erosion_cycle_record[i+1]) > steady_state_tolerance)
      {
  steady_state = false;
  return;
      }
    }
  }
  // this next one just checks if the surface topography does not change in time
  else if (steady_state_limit < 0 || current_time < steady_state_limit)
  {
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        if (abs(RasterData[i][j] - zeta_old[i][j]) > steady_state_tolerance)
  {
      steady_state = false;
      return;
  }
      }
    }
  }

  //cout << "Line 777, checking steady, iss: " << initial_steady_state << endl;
  // this is seperate logic: it simply sets the initial_steady_state flag to true
  // because if steady state has not been reached above, the return statements
  // mean this is only reached if it gets to steady state. It then checks if
  // this is the initial steady, state, and if not switches the initial_steady_state
  // flag to true
  if (initial_steady_state == false)
  {
    initial_steady_state = true;
    time_delay = current_time;

    // This is a mode of end times. 0 == default, run until the endTime
    // 1 == run until a specified time after steady state
    // 2 == run a number of cycles
    // 3 == run to steady state, then run some cycles.
    if (endTime_mode == 1 || endTime_mode == 3)
    {
      endTime += time_delay;
    }
    if (quiet == false)
    {
      cout << "\t\t\t> Initial steady state reached at " << current_time;
    }
  }
  //cout << "check iSS, it is now: " << initial_steady_state << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function checks to see if the model should record results
// If initial steady state has not been reached, recording is set to
// false: that is, the  model does not record information on the build up
// to steady state.
// JAJ, sometime 2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::check_recording( void )
{
  int num_cycles = (current_time - time_delay) / periodicity;
  if (recording)
  {
    return;
  }
  else if (initial_steady_state == false)
  {
    // If we haven't reached steady state yet, don't record any data
    recording = false;
  }
  else if (K_mode == 0 && D_mode == 0)
  {
    // If we aren't changing any of the forcing parameters, we can record
    // as soon as we hit steady state
    recording = true;
  }
  else if (num_cycles >= 1)
  {
    // If we are changing forcing parameters, we need to wait until one cycle has
    // passed, as there is a small adjustment period associated with the first cycle
    recording = true;
  }
  else
  {
    recording = false;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This checks on the ending condition of the model run
// endTime_mode:
//  1 == The end time is just some fixed time after initial steady state
//  2 == The end time is after a fixed number of cycles
//  3 == End time is after steady state, but waits for a fixed number of cycles
//       before ending
//
// Returns a boolean that is true if the end time has been reached
// and false if end time has not been reached
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDRasterModel::check_end_condition( void )
{

  //if (not quiet)
  //{
  //  cout << "LINE 888 End time mode is: " << endTime_mode
  //      << " and current time is: " << current_time << endl;
  //}
  int num_cycles;

  // if this is not cyclic, first just get a placeholder for the number of cycles
  if (K_mode != 0 || D_mode != 0)
  {
    num_cycles = cycle_number-1;
  }
  else     // if not calculate time based on current time
  {
    num_cycles = (current_time - time_delay) / periodicity;
  }
  float endTime_adjusted;

  // now check to see if end time is reached in a number of conditions
  switch (endTime_mode)
  {
    case 1:    // time specified is after reaching steady state
    {
      // end is only true if the time exceeds or is equal to the end time
      if (initial_steady_state == false || current_time <= endTime+timeStep)
        return false;
      else
        return true;
      break;
    }
    case 2:    // Number specified is a number of cycles of periodicity
    {
      if (initial_steady_state ==  false || num_cycles <= endTime)
        return false;
      else
        return true;
      break;
    }
    case 3:    // Time specified is after reaching steady state, but waits till a roughly exact number of cycles of periodicty have passed
    {
      if (ceil((endTime-time_delay)/periodicity) == 1)
      {
        endTime_adjusted = (1+ceil((endTime-time_delay) / periodicity)) * periodicity + time_delay;
      }
      else
      {
        endTime_adjusted = (ceil((endTime-time_delay) / periodicity)) * periodicity + time_delay;
      }
      //  if (not quiet)
      //    cout << "\n" << endTime << " " << time_delay << " " << periodicity << " " <<endTime_adjusted << " hi " << endl;
      endTime = endTime_adjusted;
      if (initial_steady_state == false || current_time < endTime_adjusted+timeStep)
      {
        return false;
      }
      else
      {
        return true;
      }
      break;
    }
    default:
      if (current_time >= endTime)
      {
        return true;
      }
      else
      {
        return false;
      }
      break;
  };
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function checks to see if this is a periodic run. If it is,
// it sets the times to align with the period
// JAJ, sometime 2014
// Note this only realy comes into play if period mode == 2 or 4
// period_mode means
// 1 (default) one periodicity used without
// 2 Two periodicities that switch at a given interval
// 3 Two periodicities used as a compound sin wave
// 4 Same as three, but weightings switch at a given interval (as in 2)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::check_periodicity_switch( void )
{
  // don't do anything if not periodic
  if ((K_mode == 0 && D_mode == 0) || (initial_steady_state == false && cycle_steady_check == false))
    return;
  else if (period_mode == 2 || period_mode == 4)
  {
    // this syncs the periods if there are multiple periodicities
   //if (not quiet)
    //{
    //  cout << "Periodic run. ::check_periodicity_switch. Adjusting times" << endl;
    //}
    float p = periodicity;

    float swap;
    float t;
    if (endTime_mode == 2)
      t = switch_time * p;
    else if (endTime_mode == 3)
      t = ceil((switch_time)/p) * p;
    else
      t = switch_time;

    if (current_time-time_delay > t + switch_delay)
    {
      // Time to switch periodicities yo
      /// Possible problem here, if running sequential models
      /// we can't remeber which periodicity was the original one
      swap = periodicity;
      periodicity = periodicity_2;
      periodicity_2 = swap;
      // Bump up the time til the next switch
      switch_delay = current_time - time_delay - timeStep;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// If the periodic model cycles over 100 times this returns true
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDRasterModel::check_if_hung( void )
{
  int num_cycles = (current_time) / periodicity;
  return false;
  switch (endTime_mode){
    case 1:
      if (initial_steady_state == false && current_time > endTime*100)
        return true;
      else
        return false;
      break;
    case 2:
      if (initial_steady_state == false && num_cycles > endTime * 100)
        return true;
      else
        return false;
      break;
    case 3:
      if (initial_steady_state == false && current_time > endTime*100)
        return true;
      else
        return false;
      break;
    default:
      return false;
      break;
  };
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Resets the data members total_erosion and total_response to 0
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::reset_model( void )
{
  total_erosion = 0;
  total_response = 0;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// @@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@
// Deal with the BOUNDARY CONDITIONS
// @@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// BUFFER SURFACE
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// 2 functions to buffer the raster surface.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This first function is used as a simple way to implement boundary conditions,
// particularly no flux and periodic boundary conditions.
// The buffered surface has NRows+2 rows and NCols+2 columns.
// The integer b_type sets the type of boundary conditions, but currently there
// is only one implementation: no flux across N and S; periodic for E and W.
//----------------------------------------------------------------------------
LSDRasterModel LSDRasterModel::create_buffered_surf(int b_type)
{
  Array2D<float> surf = RasterData.copy();
  Array2D<float> buff(NRows+2,NCols+2);
  Array2D<float> buff_surf = buff.copy();
//
//   switch(b_type)
//   {
//     case 1:
      // first set up the corners: these data points have no impact on calculations
      // but are here if there is some unforseen instability
      buff_surf[0][0] = surf[0][0];
      buff_surf[NRows+1][0] = surf[NRows-1][0];
      buff_surf[0][NCols+1] = surf[0][NCols-1];
      buff_surf[NRows+1][NCols+1] = surf[NRows-1][NCols-1];
      // now get the periodic boundaries
      for (int row = 0; row<NRows; row++)
      {
        buff_surf[row+1][0] = surf[row][NCols-1];
        buff_surf[row+1][NCols+1] = surf[row][0];
      }
      // now get the no flux boundaries
      for (int col = 0; col<NCols; col++)
      {
        buff_surf[0][col+1] = surf[0][col];
        buff_surf[NRows+1][col+1] = surf[NRows-1][col];
      }
      // now copy the interior
      for (int row=0; row<NRows; row++)
      {
        for (int col =0; col<NCols; col++)
        {
          buff_surf[row+1][col+1]=surf[row][col];
        }
      }
      LSDRasterModel BufferedSurface(NRows+2, NCols+2, XMinimum-DataResolution, YMinimum-DataResolution, DataResolution, NoDataValue, buff_surf, GeoReferencingStrings);
      return BufferedSurface;
//       break;
//
//     default:
//       cout << "You chose and invalid boundary condition" << endl;
//       exit(1);
//   }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//------------------------------------------------------------------------------
// This second version has periodic boundaries at E and W boundaries, and
// Neumann boundary conditions (prescribed elevations) at the N and S
// boundaries.
//------------------------------------------------------------------------------
LSDRasterModel LSDRasterModel::create_buffered_surf(float South_boundary_elevation,float North_boundary_elevation)
{
  Array2D<float> surf = RasterData.copy();
  Array2D<float> buff(NRows+2,NCols+2);
  Array2D<float> buff_surf = buff.copy();

  // now get the periodic boundaries
  for (int row = 0; row<NRows; row++)
  {
    buff_surf[row+1][0] = surf[row][NCols-1];
    buff_surf[row+1][NCols+1] = surf[row][0];
  }
  // now get the fixed elevation boundaries
  for (int col = 0; col<NCols+2; col++)
  {
    buff_surf[0][col] = South_boundary_elevation;
    buff_surf[NRows+1][col] = North_boundary_elevation;
  }
  // now copy the interior
  for (int row=0; row<NRows; row++)
  {
    for (int col =0; col<NCols; col++)
    {
      buff_surf[row+1][col+1]=surf[row][col];
    }
  }
  LSDRasterModel BufferedSurface(NRows+2, NCols+2, (XMinimum-DataResolution), (YMinimum-DataResolution), DataResolution, NoDataValue, buff_surf, GeoReferencingStrings);
  return BufferedSurface;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function seems to just set some flags within the data members
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::interpret_boundary(short &dimension, bool &periodic, int &size)
{
  dimension = 0;
  for (int i=0; i<4; ++i)
  {
    if (boundary_conditions[i][0] == 'b')
      dimension = i % 2;
  }
  if (boundary_conditions[1-dimension][0] == 'p' || boundary_conditions[3-dimension][0] == 'p')
  {
    periodic = true;
    if ( (boundary_conditions[1-dimension][0] && boundary_conditions[3-dimension][0] == 'p') == false)
      if ( quiet == false) cout << "Warning! Entered one boundary as periodic, but not t'other! Assuming both are periodic." << endl;
  }
  if (dimension == 0)
    size = (NRows-2)*NCols;
  else
    size = NRows*(NCols-2);

  if (dimension != 0 && dimension != 1)
  {
    cerr << "Warning line " << __LINE__ << ": Variable 'dimension' should have a value of 0 or 1" << endl;
    exit(1);
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function checks on nodes to see if they are base level nodes
// This is generally used in the FASTSCAPE algorithm as it starts from base
// level nodes and works its way up
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDRasterModel::is_base_level( int i, int j )
{
  // note: this only checks for base level nodes along the edges of the raster
  if (i == 0 && boundary_conditions[0][0] == 'b')
    return true;
  else if (j == 0 && boundary_conditions[3][0] == 'b')
    return true;
  else if (i == NRows-1 && boundary_conditions[2][0] == 'b')
    return true;
  else if (j == NCols-1 && boundary_conditions[1][0] == 'b')
    return true;
  else
    return false;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks on a boundary to find the maximum elevation at that boundary
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::find_max_boundary(int boundary_number)
{
  float max_val = 0;
  int i, j;
  switch (boundary_number % 2)
  {
    case 0:
      if (boundary_number == 0)
        i = 0;
      else
        i = NRows - 1;
      for (j=0; j<NCols; ++j)
      {
        if (RasterData[i][j] > max_val)
          max_val = RasterData[i][j];
      }
      break;

    case 1:
      if (boundary_number == 1)
        j = NCols - 1;
      else
        j = 0;

      for (i=0; i<NRows; ++i)
      {
        if (RasterData[i][j] > max_val)
          max_val = RasterData[i][j];
      }
      break;
  }
  return max_val;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



////------------------------------------------------------------------------------
//// impose_channels: this imposes channels onto the landscape
//// You need to print a channel to csv and then load the data
////------------------------------------------------------------------------------
void LSDRasterModel::impose_channels(LSDSpatialCSVReader& source_points_data)
{

  string column_name = "elevation(m)";


  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // need to fill the raster to ensure there are no internal base level nodes
  cout << "I am going to fill" << endl;
  float slope_for_fill = 0.0001; 
  cout << "Filling." << endl;
  LSDRaster filled_topography = temp.fill(slope_for_fill);

  cout << "Getting the flow info. This might take some time." << endl;
  LSDFlowInfo flow(boundary_conditions, filled_topography);
  // update the raster
  zeta = filled_topography.get_RasterData();

  // Get the local node index as well as the elevations
  vector<int> ni = source_points_data.get_nodeindices_from_lat_long(flow);
  vector<float> elev = source_points_data.data_column_to_float(column_name);
  // make the map
  cout << "I am making an elevation map. This will not work if points in the raster lie outside of the csv channel points." << endl;
  int row,col;
  for(int i = 0; i< int(ni.size()); i++)
  {
    flow.retrieve_current_row_and_col( ni[i], row, col);
    zeta[row][col] = elev[i];
  }


  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}


////------------------------------------------------------------------------------
//// impose_channels: this imposes channels onto the landscape
//// You need to print a channel to csv and then load the data
////------------------------------------------------------------------------------
LSDSpatialCSVReader LSDRasterModel::get_channels_for_burning(int contributing_pixels)
{

  string column_name = "elevation(m)";

  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // need to fill the raster to ensure there are no internal base level nodes
  cout << "I am going to fill" << endl;
  float slope_for_fill = 0.0001; 
  cout << "Filling." << endl;
  LSDRaster filled_topography = temp.fill(slope_for_fill);

  cout << "Getting the flow info. This might take some time." << endl;
  LSDFlowInfo FlowInfo(boundary_conditions, filled_topography);

  // calculate the flow accumulation
  cout << "\t Calculating flow accumulation (in pixels)..." << endl;
  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  //get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(FlowAcc, contributing_pixels);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  
  // print the network
  string chan_fname = "./temp_channels";
  string full_chan_fname = "./temp_channels.csv";
  ChanNetwork.PrintChannelNetworkToCSV_WithElevation(FlowInfo, chan_fname,filled_topography);

  // Now load a csv object
  LSDRasterInfo RI(temp);
  LSDSpatialCSVReader source_points_data( RI, full_chan_fname );

  return source_points_data;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// CALCULATE EROSION RATES
// Simple function that creates an array with the erosion rates for a given
// timestep, calculated by diffencing elevation rasters with consecutive
// timesteps.
//------------------------------------------------------------------------------
Array2D<float> LSDRasterModel::calculate_erosion_rates( void )
{
  // create the erosion array
  Array2D<float> ErosionRateArray(NRows,NCols,NoDataValue);

  // first check to see if zeta_old exists
  if (zeta_old.dim1() != NRows || zeta_old.dim2() != NCols)
  {
    cout << "LSDRasterModel::calculate_erosion_rates, WARNING zeta_old doesn't exist" << endl;
  }
  else
  {
    // loop through all the raster data getting erosion rate using the
    // get_erosion_at_cell data member
    for(int row=0; row<NRows; ++row)
    {
      for(int col=0; col<NCols; ++col)
      {
        if(RasterData[row][col]!=NoDataValue)
        {
          ErosionRateArray[row][col] = get_erosion_at_cell(row, col);
        }
      }
    }
  }
  return ErosionRateArray;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function does the actual calculation of the erosion rates cell by cell
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_erosion_at_cell(int i, int j)
{
  // note, this does not check to see if zeta_old exists
  return (zeta_old[i][j]-RasterData[i][j]+get_uplift_at_cell(i,j))/timeStep;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the total erosion rate over the last timestep
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_total_erosion_rate_over_timestep()
{

  float erate_total = 0;
  int N_erate = 0;
  // first check to see if zeta_old exists
  if (zeta_old.dim1() != NRows || zeta_old.dim2() != NCols)
  {
    cout << "LSDRasterModel::calculate_erosion_rates, WARNING zeta_old doesn't exist" << endl;
  }
  else
  {
    // loop through all the raster data getting erosion rate using the
    // get_erosion_at_cell data member
    for(int row=0; row<NRows; ++row)
    {
      for(int col=0; col<NCols; ++col)
      {
        if(RasterData[row][col]!=NoDataValue)
        {
          // make sure this is not a base level node
          if(is_base_level(row,col) == false)
          {
            erate_total+= get_erosion_at_cell(row, col);
            N_erate++;
          }
        }
      }
    }
  }
  return erate_total/float(N_erate);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// UPLIFT SURFACE
// Apply uplift field to the raster.  Overloaded function so that the first
// simply considers uniform uplift, the second allows user to use a prescribed
// uplift fields of greater complexity
//------------------------------------------------------------------------------
// Uniform uplift
//------------------------------------------------------------------------------
LSDRasterModel LSDRasterModel::uplift_surface(float UpliftRate, float dt)
{
  // copy the underlying surface
  Array2D<float> ZetaRaster;
  ZetaRaster = RasterData.copy();
  for(int row=0; row<NRows; ++row)
  {
    for(int col=0; col<NCols; ++col)
    {
      if(get_data_element(row,col)!=NoDataValue)
      {
        ZetaRaster[row][col] += UpliftRate*dt;
      }
    }
  }

  // make a new rastermodel with the updated data.
  // SMM Note: this is a bit risky since only a very limited subset of the
  // data members are translated across
  LSDRasterModel Zeta(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ZetaRaster, GeoReferencingStrings);
  return Zeta;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//------------------------------------------------------------------------------
// Specified uplift field
// uplift field should be specified as an array with the same dimensions as the
// elevation raster, permitting non-uniform uplift fields to be applied in the
// model.
//------------------------------------------------------------------------------
LSDRasterModel LSDRasterModel::uplift_surface(Array2D<float> UpliftRate, float dt)
{
  Array2D<float> ZetaRaster;
  ZetaRaster = RasterData.copy();
  for(int row=0; row<NRows; ++row)
  {
    for(int col=0; col<NCols; ++col)
    {
      if(get_data_element(row,col)!=NoDataValue)
      {
        ZetaRaster[row][col] += UpliftRate[row][col]*dt;
      }
    }
  }

  // make a new rastermodel with the updated data.
  // SMM Note: this is a bit risky since only a very limited subset of the
  // data members are translated across
  LSDRasterModel Zeta(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ZetaRaster, GeoReferencingStrings);
  return Zeta;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This uplift tool uses the get_uplift_at_cell member function to modify
// directly the elevation raster.
// Base level nodes are not uplifted
// Author JAJ
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::uplift_surface( void )
{
  for(int row=0; row<NRows; ++row)
  {
    for(int col=0; col<NCols; ++col)
    {
      // if it is a base level node, don't uplift
      if (is_base_level(row, col))
        continue;
      if(get_data_element(row,col)!=NoDataValue)
      {
        RasterData[row][col] += get_uplift_at_cell(row, col);
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This just gets the maximum uplift from the data members
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_max_uplift( void )
{
  return max_uplift;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This returns an uplift field in the form of an array
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D <float> LSDRasterModel::generate_uplift_field( int mode, float maximum_uplift )
{
  // Generates an uplift field from some default functions
  // Maybe worth splitting into different methods?
  uplift_mode = mode;
  max_uplift = maximum_uplift;

  Array2D <float> uplift(NRows, NCols);

  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      uplift[i][j] = get_uplift_at_cell(i, j);
    }
  }
  return uplift;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This returns an uplift field in the form of an array
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D <float> LSDRasterModel::generate_uplift_field( void )
{
  // Generates an uplift field from some default functions
  Array2D <float> uplift(NRows, NCols,0.0);

  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      uplift[i][j] = get_uplift_at_cell(i, j);
    }
  }
  return uplift;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets uplift (in distance, not a rate) at a particular cell
// It is based on the uplift mode,
// which is
// default == block uplift
// 1 == tilt block
// 2 == gaussian
// 3 == quadratic
// 4 == periodic
// The function is called by the generate_uplift_field  and uplift_surface
// member functions
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_uplift_at_cell(int i, int j)
{
  float result;

  // don't do anything if this is a base level node
  if (is_base_level(i,j))
    return 0;
  else
    switch (uplift_mode)
    {
      case 1:    // Tilt block (SMM: seems to only tilt in one direction)
      {
        if (baseline_uplift > 0 && baseline_uplift < max_uplift)
        {
          // this tilt block is maximum in row 1 and at the baseline_uplift
          // in row NRows-2
          float m = (baseline_uplift-get_max_uplift())/(NRows-3);
          result = float(i)*m+get_max_uplift()-m;
        }
        else
        {
          result = (NRows - i - 1) * get_max_uplift() / ((float) NRows - 1);
        }
        break;
      }
      case 2:    // Gausian
      {
        // Gaussian parameters
        int mu_i = NRows/2;
        int mu_j = NCols/2;
        float sigma_i = NRows/10;
        float sigma_j = NCols/10;

        result = get_max_uplift()*pow(1.1, -((i-mu_i)*(i-mu_i)/(2*sigma_i*sigma_i) + (j-mu_j)*(j-mu_j)/(2*sigma_j*sigma_j) ));
        break;
      }
      case 3:    // polynomial
      {
        result = get_max_uplift() * ( -pow((2.0*i/(NRows-1) - 1),2) - pow((2.0*j/(NCols-1) - 1), 2) + 1);
        if (result < 0)
          result = 0;
        break;
      }
      case 4:
      {
        result = periodic_parameter( get_max_uplift(), uplift_amplitude );
        break;
      }
      default:
        result = get_max_uplift();
        break;
    }
  // the uplift is multiplied by the timestep so uplift is passed as a distance
  // and not a rate.
  return result * timeStep;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets uplift rate at a particular cell
// It is based on the uplift mode,
// which is
// default == block uplift
// 1 == tilt block
// 2 == gaussian
// 3 == quadratic
// 4 == periodic
// The function is called by the generate_uplift_field  and uplift_surface
// member functions
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_uplift_rate_at_cell(int i, int j)
{
  float result;


  // don't do anything if this is a base level node
  if (is_base_level(i,j))
    return 0;
  else
    switch (uplift_mode)
    {
      case 1:  // Tilt block (SMM: seems to only tilt in one direction)
      {
        if (baseline_uplift > 0 && baseline_uplift < max_uplift)
        {
          // this tilt block is maximum in row 1 and at the baseline_uplift
          // in row NRows-2
          float m = (baseline_uplift-get_max_uplift())/(NRows-3);
          result = float(i)*m+get_max_uplift()-m;
        }
        else
        {
          result = (NRows - i - 1) * get_max_uplift() / ((float) NRows - 1);
        }
        break;
      }
      case 2:    // Gausian
      {
        // Gaussian parameters
        int mu_i = NRows/2;
        int mu_j = NCols/2;
        float sigma_i = NRows/10;
        float sigma_j = NCols/10;

        result = get_max_uplift()*pow(1.1, -((i-mu_i)*(i-mu_i)/(2*sigma_i*sigma_i) + (j-mu_j)*(j-mu_j)/(2*sigma_j*sigma_j) ));
        break;
      }
      case 3:    // polynomial
      {
        result = get_max_uplift() * ( -pow((2.0*i/(NRows-1) - 1),2) - pow((2.0*j/(NCols-1) - 1), 2) + 1);
        if (result < 0)
          result = 0;
        break;
      }
      case 4:
      {
        result = periodic_parameter( get_max_uplift(), uplift_amplitude );
        break;
      }
      case 5:
      {
        // The uplift field is a distance, so to get the rate we need to divide by the timestep.
        // This is inherited from an old way of calculating uplift.
        // I will need to fix in future but don't know what all the dependencies are
        // so can't do it now.
        result = uplift_field[i][j]/timeStep;

      }
      default:
        result = get_max_uplift();
        break;
    }
  // the uplift is multiplied by the timestep so uplift is passed as a distance
  // and not a rate.
  return result;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Checks to see if the uplift field is the correct size
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::check_and_correct_uplift_field()
{
  cout << "I am checking to make sure the uplift field is the same dimensions" << endl;
  cout << " as the model data. " << endl;

  cout << "NRows: " << NRows << " NCols: " << NCols << endl;
  cout << "Data Rows: " << RasterData.dim1() << " NCols: " << RasterData.dim2() << endl;
  cout << "Uplift rows: " << uplift_field.dim1() << " uplift cols: " << uplift_field.dim2() << endl;

  // check to see, if the uplift mode is based on an uplift field,
  // that the uplift field is the same size as the raster
  if(uplift_mode == 5)
  {
    cout << "You have decided to uyse the uplift field matrix! " << endl;
  }


}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the average uplift rate. It excludes N and S boundaries.
// Can be used to get tectonic fluxes
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_average_upflit_rate_last_timestep()
{

  float avg_uplift_rate;
  float tot_urate = 0;
  int N_U = 0;

  for (int row = 0; row<NRows-1; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if(RasterData[row][col] != NoDataValue)
      {
        // check to see if it is a base level cell
        if(is_base_level(row,col) == false)
        {
          tot_urate += get_uplift_rate_at_cell(row, col);
          N_U++;
        }
      }
    }
  }

  avg_uplift_rate = tot_urate/float(N_U);
  return avg_uplift_rate;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// set_uplift_field_to_block_uplift
// This function sets the uplift_field data member to block uplift at a given
// uplift rate
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::set_uplift_field_to_block_uplift(float uplift_rate)
{
  Array2D<float> uplift_rate_field(NRows,NCols,uplift_rate);
  uplift_field = uplift_rate_field.copy();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// CREATE PRECIPITION FLUX ARRAY
// Produces precipitation array from provided precipitation rate.
//------------------------------------------------------------------------------
Array2D<float> LSDRasterModel::precip_array_from_precip_rate(float precip_rate)
{
  float precip_flux = DataResolution*DataResolution*precip_rate;
  Array2D<float> precip_start(NRows,NCols,precip_flux);
  Array2D<float> PrecipitationFluxArray = precip_start.copy();
  return PrecipitationFluxArray;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// TOPOGRAPHIC DERIVATIVES
//----------------------------------------------------------------------------
// get_slopes
//----------------------------------------------------------------------------
// Specifically, this function gets the topographic slopes, as required for
// the sediment flux calculations.  The slopes are stored as two matrices, one
// that stores slopes between rows, the other which for slopes between
// columns.  Note that this is a finite volume model that utilises cubic model
// voxels. Sediment fluxes are only permitted through the faces.
//
// For slopes between columns, the entry at S[row][col] refers to the slope
// between zeta at node [row][col] and at node [row][col+1].  Likewise for the
// slopes between rows.  In short, the center points of the slopes are offset
// by 1/2 a node spacing in the positive direction.
// Note that there are NCols +1 and NRows +1 columns and rows respectively
//----------------------------------------------------------------------------
void LSDRasterModel::get_slopes(Array2D<float>& SlopesBetweenRows, Array2D<float>& SlopesBetweenCols)
{
  Array2D<float> buff_zeta = RasterData.copy();
  float inv_dx = 1/DataResolution;
  float inv_dy = 1/DataResolution;
  //cout << "LINE 50 flux_funcs, n_rows: " << n_rows
  //     << " and n_cols: " << n_cols << endl;
  //cout << "LINE 52 flux_funcs, zeta: " << zeta << endl;

  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<=NCols; col++)
    {
      //cout << "row: " << row << " and col: " << col << endl;
      SlopesBetweenCols[row][col] = (buff_zeta[row+1][col+1] - buff_zeta[row+1][col])*inv_dx;
    }
  }

  for (int row = 0; row<=NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      SlopesBetweenRows[row][col] = (buff_zeta[row+1][col+1] - buff_zeta[row][col+1])*inv_dy;
    }
  }
}

//----------------------------------------------------------------------------
// get_topographic_divergence
// gets the topographic divergence at each point in the model domain.  Use
// buffered topography
//----------------------------------------------------------------------------
Array2D<float> LSDRasterModel::get_topographic_divergence()
{
  Array2D<float> buffered_topo = RasterData.copy();
  Array2D<float> empty_div(NRows,NCols,0.0);
  Array2D<float> div_zeta = empty_div.copy();
  float s1,s2;
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      s1 = (buffered_topo[row+1][col+2]-buffered_topo[row+1][col])*0.5/DataResolution;
      s2 = (buffered_topo[row+2][col+1]-buffered_topo[row][col+1])*0.5/DataResolution;
      div_zeta[row][col] = sqrt(s1*s1+s2*s2);
    }
  }
  Array2D<float> TopoDivergence = div_zeta.copy();
  return TopoDivergence;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// HYDROLOGICAL TOOLS
//----------------------------------------------------------------------------
// calculate_channel_width_wolman
// This function calculates channel width using the wolman method.
// NOTE: typically Q_w will be in m^3/s.
// EXAMPLE: in Salmon River, Idaho (Emmett, 1975 cited in Knighton 1988):
//          k_w = 2.77 and b = 0.56. b is often assumed to be 0.5
//----------------------------------------------------------------------------
float LSDRasterModel::calculate_channel_width_wolman(float Q_w, float k_w, float b)
{
  float ChannelWidth;
  if (b == 1.0)
  {
    ChannelWidth = Q_w*k_w;
  }
  else if (b == 0.5)
  {
    ChannelWidth = k_w*sqrt(Q_w);
  }
  else
  {
    ChannelWidth = k_w*pow(Q_w,b);
  }
  return ChannelWidth;
}
//----------------------------------------------------------------------------
// array_channel_width_wolman
// this function calcualtes channel width in a stand alone module so the widths
// can be tested
//----------------------------------------------------------------------------
Array2D<float> LSDRasterModel::array_channel_width_wolman(Array2D<float>& Q_w, float& k_w, float& b)
{
  // reset the channel width array
  Array2D<float> empty_w(NRows,NCols,1.0);
  Array2D<float> ChannelWidth = empty_w.copy();

  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      ChannelWidth[row][col] = calculate_channel_width_wolman(Q_w[row][col], k_w, b);
    }
  }
  return ChannelWidth;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// EROSION RATES/SEDIMENT FLUXES
//
//----------------------------------------------------------------------------
// this caluclates the fluvial erosion rate at each point
//----------------------------------------------------------------------------
Array2D<float> LSDRasterModel::calculate_fluvial_erosion_rate(Array2D<float> ChannelWidth, Array2D<float> Q_w,
              Array2D<float> TopoDivergence, float K, float n, float m, float eros_thresh)
{
  // set up an array to populate with erosion
  Array2D<float> empty_fluv_eros(NRows,NCols,0.0);
  Array2D<float> FluvialErosionRate = empty_fluv_eros.copy();

  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      FluvialErosionRate[row][col] = K*(ChannelWidth[row][col]/DataResolution)*pow(TopoDivergence[row][col],n)*pow(Q_w[row][col],m) - eros_thresh;
      if (FluvialErosionRate[row][col] < 0)
      {
        FluvialErosionRate[row][col] = 0;
      }
    }
  }
  return FluvialErosionRate;
}
//------------------------------------------------------------------------------

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// IMPLICIT MODEL COMPONENTS
//------------------------------------------------------------------------------
// Implicit schemes for combination of hillslope sediment transport using
// non-linear hillslope transport law, and fluvial erosion.  This is essentially
// the implicit implementation of MuddPILE, but has been modified so that now
// fluvial erosion is undertaken using FASTSCAPE (Braun and Willet, 2013), which
// greatly increases computational efficiency.
//------------------------------------------------------------------------------
// calculate_k_values_for_assembly_matrix/mtl_initiate_assembler_matrix
// this function creates vectors of integers that refer to the k values, that is
// the index into the vectorized matrix of zeta values, that is used in the assembly matrix
// the number of elements in the k vectors is N_rows*N_cols
//------------------------------------------------------------------------------
void LSDRasterModel::calculate_k_values_for_assembly_matrix(int NRows, int NCols, vector<int>& k_value_i_j,
                      vector<int>& k_value_ip1_j,  vector<int>& k_value_im1_j, vector<int>& k_value_i_jp1,
                      vector<int>& k_value_i_jm1)
{
  int N_elements_in_k_vec = NRows*NCols;

  // initialize the vectors with empty values
  vector<int> empty_vec(N_elements_in_k_vec,0);
  k_value_i_j   = empty_vec;
  k_value_ip1_j = empty_vec;
  k_value_im1_j = empty_vec;
  k_value_i_jp1 = empty_vec;
  k_value_i_jm1 = empty_vec;

  // we loop through each node
  int counter = 0;
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      k_value_ip1_j[counter] = NCols*(row+2)+col;
      k_value_im1_j[counter] = NCols*row+col;
      k_value_i_j[counter] = NCols*(row+1)+col;

      // logic for west periodic boundary
      if(col == 0)
      {
        k_value_i_jp1[counter] = NCols*(row+1)+col+1;
        k_value_i_jm1[counter] = NCols*(row+1)+NCols-1;
      }
      // logic for east periodic boundary
      else if(col == NCols-1)
      {
        k_value_i_jp1[counter] = NCols*(row+1);
        k_value_i_jm1[counter] = NCols*(row+1)+col-1;

      }
      // logic for rest of matrix
      else
      {
        k_value_i_jp1[counter] = NCols*(row+1)+col+1;
        k_value_i_jm1[counter] = NCols*(row+1)+col-1;
      }

      // increment counter
      counter++;
    }
  }
}

void LSDRasterModel::mtl_initiate_assembler_matrix(int& problem_dimension,
             float& inv_dx_S_c_squared, float& inv_dy_S_c_squared, float& dx_front_term,
               float& dy_front_term, vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,
       vector<int>& vec_k_value_im1_j, vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1)
{
  float dx = DataResolution;
  float dy = DataResolution;
  float D = get_D();

  inv_dx_S_c_squared = 1/(dx*dx*S_c*S_c);
  inv_dy_S_c_squared = 1/(dy*dy*S_c*S_c);
  dx_front_term = timeStep*D/(dx*dx);
  dy_front_term = timeStep*D/(dy*dy);

  problem_dimension = (NRows+2)*NCols;
  calculate_k_values_for_assembly_matrix(NRows, NCols, vec_k_value_i_j, vec_k_value_ip1_j,
                      vec_k_value_im1_j, vec_k_value_i_jp1, vec_k_value_i_jm1);

}

//------------------------------------------------------------------------------
// mtl_assemble_matrix
// this function assembles the solution matrix for nonlinear creep transport
//------------------------------------------------------------------------------
void LSDRasterModel::mtl_assemble_matrix(Array2D<float>& zeta_last_iter, Array2D<float>& zeta_last_timestep,
             Array2D<float>& zeta_this_iter, Array2D<float>& uplift_rate, Array2D<float>& fluvial_erosion_rate,
             mtl::compressed2D<float>& mtl_Assembly_matrix, mtl::dense_vector<float>& mtl_b_vector,
             float dt, int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
             float dx_front_term, float dy_front_term,
             float South_boundary_elevation, float North_boundary_elevation,
             vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,vector<int>& vec_k_value_im1_j,
             vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1)
{
  // the coefficients in the assembly matrix
  float A,B,C,D;

  // reset the assembly and b vector
  mtl_Assembly_matrix = 0.0;
  mtl_b_vector = 0.0;

  // create the inserter. This is deleted when this function is exited
  mtl::mat::inserter< mtl::compressed2D<float> > ins(mtl_Assembly_matrix);

  // first we assemble the boundary nodes. First the nodes in row 0 (the south boundary)
  for (int k = 0; k<NCols; k++)
  {
    ins[k][k] << 1.0;
    mtl_b_vector[k] =  South_boundary_elevation;//zeta_last_timestep[0][k];
  }

  // now assemble the north boundary
  int starting_north_boundary = (NRows+1)*(NCols);
  int one_past_last_north_boundary = (NRows+2)*NCols;
  for (int k = starting_north_boundary; k < one_past_last_north_boundary; k++)
  {
    ins[k][k] << 1.0;
    mtl_b_vector[k] = North_boundary_elevation;//zeta_last_iter[NCols-1][k];
  }

  // create the zeta matrix that includes the boundary conditions
  Array2D<float> zeta_for_implicit(NRows+2, NCols,0.0);
  for (int col = 0; col<NCols; col++)
  {
    zeta_for_implicit[0][col] = zeta_last_iter[0][col];
    zeta_for_implicit[NRows+1][col] = zeta_last_iter[NRows-1][col];
  }
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      zeta_for_implicit[row+1][col] = zeta_last_iter[row][col];
    }
  }

  // now assemble the rest
  // we loop through each node
  int counter = 0;
  float b_value;
  int k_value_i_j,k_value_ip1_j,k_value_im1_j,k_value_i_jp1,k_value_i_jm1;
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if  (col == 0 || col == NCols-1)
      {
        b_value = zeta_last_iter[row][col];
      }
      else
      {
        b_value = zeta_last_timestep[row][col]+dt*uplift_rate[row][col]-dt*fluvial_erosion_rate[row][col];
      }
      k_value_ip1_j = vec_k_value_ip1_j[counter];
      k_value_im1_j = vec_k_value_im1_j[counter];
      k_value_i_j   = vec_k_value_i_j[counter];
      k_value_i_jp1 = vec_k_value_i_jp1[counter];
      k_value_i_jm1 = vec_k_value_i_jm1[counter];

      A =  dy_front_term/(1 -
              (zeta_for_implicit[row+2][col]-zeta_for_implicit[row+1][col])*
              (zeta_for_implicit[row+2][col]-zeta_for_implicit[row+1][col])*
          inv_dy_S_c_squared);
      B = dy_front_term/(1 -
              (zeta_for_implicit[row+1][col]-zeta_for_implicit[row][col])*
              (zeta_for_implicit[row+1][col]-zeta_for_implicit[row][col])*
          inv_dy_S_c_squared);

      // logic for west periodic boundary
      if(col == 0)
      {
        C = dx_front_term/(1 -
                 (zeta_for_implicit[row+1][col+1]-zeta_for_implicit[row+1][col])*
                 (zeta_for_implicit[row+1][col+1]-zeta_for_implicit[row+1][col])*
             inv_dx_S_c_squared);
        D = dx_front_term/(1 -
                 (zeta_for_implicit[row+1][col]-zeta_for_implicit[row+1][NCols-1])*
                 (zeta_for_implicit[row+1][col]-zeta_for_implicit[row+1][NCols-1])*
             inv_dx_S_c_squared);
         //A=B=C=D=0;
      }
      // logic for east periodic boundary
      else if(col == NCols-1)
      {
        C = dx_front_term/(1 -
                 (zeta_for_implicit[row+1][0]-zeta_for_implicit[row+1][col])*
                 (zeta_for_implicit[row+1][0]-zeta_for_implicit[row+1][col])*
             inv_dx_S_c_squared);
        D = dx_front_term/(1 -
                 (zeta_for_implicit[row+1][col]-zeta_for_implicit[row+1][col-1])*
                 (zeta_for_implicit[row+1][col]-zeta_for_implicit[row+1][col-1])*
             inv_dx_S_c_squared);
        //A=B=C=D=0;

      }
      // logic for rest of matrix
      else
      {
        C = dx_front_term/(1 -
                 (zeta_for_implicit[row+1][col+1]-zeta_for_implicit[row+1][col])*
                 (zeta_for_implicit[row+1][col+1]-zeta_for_implicit[row+1][col])*
             inv_dx_S_c_squared);

        D = dx_front_term/(1 -
                 (zeta_for_implicit[row+1][col]-zeta_for_implicit[row+1][col-1])*
                 (zeta_for_implicit[row+1][col]-zeta_for_implicit[row+1][col-1])*
             inv_dx_S_c_squared);

      }

      // place the values in the assembly matrix and the b vector
      mtl_b_vector[k_value_i_j] = b_value;
      ins[k_value_i_j][k_value_ip1_j] << -A;
      ins[k_value_i_j][k_value_im1_j] << -B;
      ins[k_value_i_j][k_value_i_jp1] << -C;
      ins[k_value_i_j][k_value_i_jm1] << -D;
      ins[k_value_i_j][k_value_i_j] << 1+A+B+C+D;

      counter++;
    }
  }
}

//------------------------------------------------------------------------------
// mtl_solve_assembler_matrix
// this function assembles the solution matrix
//------------------------------------------------------------------------------
void LSDRasterModel::mtl_solve_assembler_matrix(Array2D<float>& zeta_last_iter, Array2D<float>& zeta_last_timestep,
             Array2D<float>& zeta_this_iter, Array2D<float>& uplift_rate, Array2D<float>& fluvial_erosion_rate,
             float dt, int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
             float dx_front_term, float dy_front_term,
             vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j, vector<int>& vec_k_value_im1_j,
             vector<int>& vec_k_value_i_jp1, std::vector<int>& vec_k_value_i_jm1,
             float South_boundary_elevation, float North_boundary_elevation)
{
  // reset the zeta array for this iteration
  Array2D<float> empty_zeta(NRows,NCols,0.0);
  zeta_this_iter = empty_zeta.copy();
  //zeta_this_iter(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, empty_zeta);
  // create a mtl matrix
  // NOTE: you could probably save time by creating the mtl matrix and vector
  // in main()
  mtl::compressed2D<float> mtl_Assembly_matrix(problem_dimension, problem_dimension);
  mtl::dense_vector<float> mtl_b_vector(problem_dimension,0.0);

  // assemble the matrix
  mtl_assemble_matrix(zeta_last_iter, zeta_last_timestep, zeta_this_iter,
            uplift_rate, fluvial_erosion_rate, mtl_Assembly_matrix, mtl_b_vector,
            dt, problem_dimension, inv_dx_S_c_squared, inv_dy_S_c_squared,
            dx_front_term, dy_front_term, South_boundary_elevation, North_boundary_elevation,
            vec_k_value_i_j, vec_k_value_ip1_j, vec_k_value_im1_j, vec_k_value_i_jp1, vec_k_value_i_jm1);
  //std::cout << "matrix assembled!" << endl;

  //std::ofstream assembly_out;
  //assembly_out.open("assembly.data");
  //assembly_out << mtl_Assembly_matrix << endl;
  //assembly_out.close();

  // now solve the mtl system
  // Create an ILU(0) preconditioner
  bool show_time = false;
  long time_start, time_end, time_diff;
  time_start = time(NULL);
  itl::pc::ilu_0< mtl::compressed2D<float> > P(mtl_Assembly_matrix);
  mtl::dense_vector<float> mtl_zeta_solved_vector(problem_dimension);
  itl::basic_iteration<float> iter(mtl_b_vector, 500, 1.e-8);
  bicgstab(mtl_Assembly_matrix, mtl_zeta_solved_vector, mtl_b_vector, P, iter);
  time_end = time(NULL);
  time_diff = time_end-time_start;
  if (show_time)
  {
    std::cout << "iter MTL bicg took: " << time_diff << endl;
  }

  // now reconstitute zeta
  int counter = 0;//NCols;
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      zeta_this_iter[row][col] = mtl_zeta_solved_vector[counter];
      counter++;
    }
  }
}

//------------------------------------------------------------------------------
// nonlinear_creep_timestep
// do a creep timestep.  This function houses the above two functions to
// undertake model timestep using implicit implementation of the nonlinear
// transport law.
// NOTE you need to run mtl_initiate_assembler_matrix before you run this function
//------------------------------------------------------------------------------
void LSDRasterModel::nonlinear_creep_timestep(Array2D<float>& fluvial_erosion_rate, float iteration_tolerance,
      int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
      float dx_front_term, float dy_front_term, vector<int>& vec_k_value_i_j,
      vector<int>& vec_k_value_ip1_j, vector<int>& vec_k_value_im1_j,
      vector<int>& vec_k_value_i_jp1,  vector<int>& vec_k_value_i_jm1,
      float South_boundary_elevation, float North_boundary_elevation)
{
   // reset zeta_old and zeta_intermediate
  Array2D<float> zeta = RasterData.copy();
  Array2D<float> zeta_old = zeta.copy();
  Array2D<float> zeta_intermediate = zeta.copy();

  // set up residual
  float residual;
  float N_nodes = float(NRows*NCols);
  int iteration = 0;
  int Max_iter = 100;
  do
  {
    residual = 0.0;
    mtl_solve_assembler_matrix(zeta, zeta_old, zeta_intermediate,
        uplift_field, fluvial_erosion_rate, timeStep,
        problem_dimension, inv_dx_S_c_squared, inv_dy_S_c_squared,
        dx_front_term, dy_front_term, vec_k_value_i_j, vec_k_value_ip1_j,
        vec_k_value_im1_j, vec_k_value_i_jp1, vec_k_value_i_jm1,
        South_boundary_elevation, North_boundary_elevation);


    // check the residuals (basically this is the aveage elevation change between intermediate
    // zeta values
    for (int row = 0; row<NRows; row++)
    {
      for (int col = 0; col<NCols; col++)
      {
        residual+= sqrt( (zeta_intermediate[row][col]-zeta[row][col])*
                 (zeta_intermediate[row][col]-zeta[row][col]) );
      }
    }
    residual = residual/N_nodes;

    // reset last zeta
    zeta = zeta_intermediate.copy();

    iteration++;

    if (iteration%5 == 0)
    {
      //std::cout << "iteration is: " << iteration << " and residual RMSE is: " << residual << endl;
    }
    if (iteration > Max_iter)
    {
      iteration_tolerance = iteration_tolerance*10;
      iteration = 0;
    }
  } while (residual > iteration_tolerance);

  //   LSDRasterModel ZetaNew(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  // Update Raster Data
  RasterData = zeta.copy();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// RUN MODEL
//------------------------------------------------------------------------------
// A series of wrapper functions that implement the numerical model
//------------------------------------------------------------------------------
// implicit_hillslope_and_fluvial
// This function sets up a landscape evolution model run incorporating fluvial
// erosion and hillslope erosion via non-linear creep.  It calls the implicit
// implementation.
// The user should provide the parameter file which sets out the details of the
// model run.
//------------------------------------------------------------------------------
LSDRasterModel LSDRasterModel::run_model_implicit_hillslope_and_fluvial(string param_file)
{
  // parameters. All of these are initialized using the .param file
  // time and printing paramters
  float dt = timeStep;        // time spacing
  float EndTime = endTime;          // the time at which the model ends
  float PrintInterval;     // the frequency at which the model prints data

  // fluvial parameters
  float k_w;                // parameter for determining channel width from discharge in m^3/s
                            // based in Emmett (1975) as cited in Knighton (1988)from Salmon River, ID
  float b;                  // channel width exponent.
  float m;                  // area exponent for fluvial incision
  float n;                  // slope exponent for fluvial incision
  float K;                  // coefficient for fluvial incision
  float ErosionThreshold;  // threshold amunt of fluvial action required to erode bed (at this stage
                            // it is an erosion rate subtracted from the main stream power erosion rate)

  // creep-like parameters
  float K_nl;        // diffusivity of hillslope sediment
  float S_c;          // critical slope

  // forcing parameters
  float uplift_rate;      // rate of rock uplift
  float precip_rate;      // rate of precipitation

  // boundary_conditions
  float North_boundary_elevation;    // the elevation of the channels at the north
  float South_boundary_elevation;    // and south bounding nodes. Note these nodes only appear
                                      // in the buffered grid
  // file names
  string run_name;              // the name of the run
              // the file formats are: file_type.run_name.time_step.data
  //string surf_fname;              // the name of the surface file
  string area_name;
  string div_name;
  string w_name;
  string erosion_rate_fname;

  // data elements: vectors and arrays
  LSDRasterModel Zeta(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, get_RasterData(), GeoReferencingStrings);
  Array2D<float> ZetaOld;            // surface from last timestep (for getting sediment flux)
  Array2D<float> ZetaTemp(NRows,NCols);
  Array2D<float> ZetaDivergence;      // del dot zeta
  Array2D<float> PrecipitationFlux;  // precipitation flux (will be in units m^3 per second
  Array2D<float> Q_w;                // discahrge (in m^3/second)
  Array2D<float> w;                  // channel width in m
  Array2D<float> FluvialErosionRate;  // the fluvial erosion rate in m/yr
  Array2D<float> SlopesBetweenCols;
  Array2D<float> SlopesBetweenRows;
  Array2D<float> ErosionRate;

  // file for printing timesteps
  ofstream ts_out;
  cout << "LINE " << __LINE__ << ": Initializing model" << endl;

  // initiate the model
  Zeta.initialize_model(param_file, run_name, dt, EndTime, PrintInterval,
    k_w, b, m, n, K, ErosionThreshold, K_nl, S_c, uplift_rate, precip_rate,
    North_boundary_elevation, South_boundary_elevation,
    PrecipitationFlux, SlopesBetweenRows, SlopesBetweenCols, ErosionRate);

  cout << "LINE " << __LINE__ << ": Model initialized" << endl;
  Array2D<float> UpliftRate(NRows,NCols,uplift_rate);

  // print the initial condition
  float t_ime = 0;
  Array2D<float> Q_w_temp(NRows,NCols,0.0);
  Q_w = Q_w_temp.copy();
  //int print_counter = 0;
//   print_at_print_interval(run_name, t_ime, dt, N_rows, N_cols, print_interval,           // NEED TO CHANGE THIS MODULE SO THAT IT SIMPLY CALLS LSDRasterModel.write_raster()
//                 print_counter, ts_out, zeta, erosion_rate, fluvial_erosion_rate, Q_w);

//  cout << "initial_condition printed " << endl;

  // now use the tolerance method
  // initialize some variables needed to speed up the calucaltions
  int problem_dimension;
  float inv_dx_S_c_squared, inv_dy_S_c_squared, dx_front_term, dy_front_term;
  vector<int> vec_k_value_i_j;
  vector<int> vec_k_value_ip1_j;
  vector<int> vec_k_value_im1_j;
  vector<int> vec_k_value_i_jp1;
  vector<int> vec_k_value_i_jm1;
  float iteration_tolerance = 0.01;

  // now initialize the assembly matrix and some constants
  Zeta.mtl_initiate_assembler_matrix(problem_dimension,
             inv_dx_S_c_squared, inv_dy_S_c_squared, dx_front_term, dy_front_term,
               vec_k_value_i_j, vec_k_value_ip1_j, vec_k_value_im1_j, vec_k_value_i_jp1, vec_k_value_i_jm1);

  if (quiet == false) cout << "LINE " << __LINE__ << ": assembler matrix initialized" << endl;

  // do a time loop
  // now do the time loop
  while (t_ime < EndTime)
  {
    t_ime+= dt;

    if (quiet == false) cout << flush << "time is: " << t_ime << "\r";

//     // buffer the landscape
//     ZetaBuff = Zeta.create_buffered_surf(South_boundary_elevation,North_boundary_elevation);
//
//     // Fill Buffered array
//     float min_slope = 0.0001;
//     ZetaBuff.fill_overwrite(min_slope);
//     Array2D<float> ZetaUpdate(NRows,NCols);
//     // now update original surface
//     for (int i = 0; i<NRows; ++ i)
//     {
//       for (int j = 0; j<NCols; ++j)
//       {
//         ZetaUpdate[i][j]=ZetaBuff.get_data_element(i+1,j+1);
//       }
//     }
//
//     Zeta.RasterData = ZetaUpdate;

    // ERODE FLUVIALLY USING FASTSCAPE
    // Do flow routing using FlowInfo to calculate Q_w

    // now get the channel widths
    //w = array_channel_width_wolman(Q_w, k_w, b);

    // calcualte the topographic divergence
    //ZetaDivergence = Zeta.get_topographic_divergence();

    // calculate fluvial erosion rate based on Fastscape
    //FluvialErosionRate = calculate_fluvial_erosion_rate(w, Q_w, ZetaDivergence, K, n, m, ErosionThreshold);
    Array2D<float> fluvial_temp(NRows,NCols,0.0);
    //FluvialErosionRate = fluvial_temp.copy;
    // update elevations to reflect fluvial incision during timestep?
    // - we are going to need to think about this if we are to maintain stable hillslopes

  // HILLSLOPE EROSION
    // do a nonlinear creep timestep
//    cout << "/n non_linear_creep_timestep" << endl;
    Zeta.nonlinear_creep_timestep(fluvial_temp,
          iteration_tolerance, problem_dimension,
          inv_dx_S_c_squared, inv_dy_S_c_squared, dx_front_term, dy_front_term,
          vec_k_value_i_j, vec_k_value_ip1_j, vec_k_value_im1_j,
          vec_k_value_i_jp1, vec_k_value_i_jm1,
            South_boundary_elevation, North_boundary_elevation);

    // calculate erosion rate
//    ErosionRate = Zeta.calculate_erosion_rates(ZetaOld,dt);

//    print_at_print_interval(run_name, t_ime, dt, NRows, NCols, print_interval,
//            print_counter, ts_out, zeta, erosion_rate, fluvial_erosion_rate, Q_w);
  }
  ts_out.close();
  //LSDRasterModel(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Zeta.get_RasterData());
  return Zeta;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is a wrapper function used to drive a model run
// it checks parameters and flags from the data members
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::run_components( void )
{
  //recording = false;
  cycle_number = 1;
  total_erosion = 0;
  max_erosion = 0;
  min_erosion = -99;
  switch_delay = 0;
  time_delay = 0;

  stringstream ss, ss_root;

  // set the fram to the current frame
  int frame = current_frame;
  int print = 1;
  do
  {
    // Check if hung
    if ( check_if_hung() )
    {
      cout << "Model took too long to reach steady state, assumed to be stuck" << endl;
      break;
    }

    // ONLY IMPORTANT IF USING TWO PERIODICITIES
    check_periodicity_switch();     // This checks to see if this is a periodic
                                    // run, if so it adjusts time and end time
                                    // it really only comes into play if you
                                    // have more than one periodicity that
                                    // between which the model oscillates.

    // Record current topography
    zeta_old = RasterData.copy();

    // run active model components
    // first diffuse the hillslopes. Currently two options. Perhaps a flag could be
    // added so that we don't have to keep adding if statements as more
    // hillslope rules are added?
    if (hillslope)
    {
      if (nonlinear)
      {
        //soil_diffusion_fv_nonlinear();
        MuddPILE_nl_soil_diffusion_nouplift();
      }
      else
      {
        soil_diffusion_fd_linear();
      }
    }

    // sediment will have moved into channels. This is assumed to
    // be immediately removed by fluvial processes, so we use the
    // 'wash out' member function
    wash_out();

    // now for fluvial erosion
    if (fluvial)
    {
      // currently the opnly option is stream power using the FASTSCAPE
      // algorithm so there are no choices here
      fluvial_incision();
    }

    // now for isostacy. As of 17/6/2014 this wasn't working
    if (isostasy)
    {
      if (flexure)
      {
        //flexural_isostasy( 0.00001 );
        flexural_isostasy_alt( );
      }
      else
      {
        Airy_isostasy();
      }
    }

    // after everything has been eroded/deposited, the surface is uplifted
    // this currently is the block uplift option
    uplift_surface();

    // this writes a file if it is on the print interval, and writes some erosion
    // data to an array always so that steady state may be tested
    write_report();

    current_time += timeStep;
    //cout << "current_time is: " << current_time << endl;

    //cout << "Print is " << print << " and print interval is: " << print_interval << endl;
    if ( print_interval > 0 && (print % print_interval) == 0)  // write at every print interval
    {
      //cout << "Printing frame " << frame << " at time " << current_time << " years" << endl;
      print_rasters( frame );
      ++frame;
    }
    if (quiet == false) cout << "\rTime: " << current_time << " years" << flush;
    ++print;

    // check to see if steady state has been achieved
    //cout << "Line 2222, checking steady, initial steady: " << initial_steady_state << endl;
    check_steady_state();
    //cout << "Line 2224, checked, iss: " << initial_steady_state << endl;

  } while (check_end_condition() == false);

  if ( print_interval == 0 || (print_interval > 0 && ((print-1) % print_interval) != 0))
  {
    // if the print interval is 0, we set the frame to zero
    if (print_interval == 0)
    {
      frame = 0;
    }
    //if (not quiet) cout << current_time << ": " << endl;
    print_rasters( frame );
  }

  // reset the current frame
  current_frame = frame;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is a wrapper function used to drive a model run
// it checks parameters and flags from the data members
//
// Similar to run_componets but instead of running fluvial
// and uplift in seperate steps it inserts the fluvial
// and uplift matrices into the nonlinear hillslope
// solver.
//
// SMM, 06/07/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::run_components_combined( void )
{
  //recording = false;
  cycle_number = 1;
  total_erosion = 0;
  max_erosion = 0;
  min_erosion = -99;
  switch_delay = 0;
  time_delay = 0;

  stringstream ss, ss_root;

  // some fields for uplift and fluvial incision
  Array2D<float> uplift_field;
  Array2D<float> fluvial_incision_rate_field;

  // set the frame to the current frame plus 1
  int frame = current_frame;
  int print = 1;
  do
  {
    // Check if hung
    if ( check_if_hung() )
    {
      cout << "Model took too long to reach steady state, assumed to be stuck" << endl;
      break;
    }

    // ONLY IMPORTANT IF USING TWO PERIODICITIES
    check_periodicity_switch();     // This checks to see if this is a periodic
                                    // run, if so it adjusts time and end time
                                    // it really only comes into play if you
                                    // have more than one periodicity that
                                    // between which the model oscillates.

    // Record current topography
    zeta_old = RasterData.copy();

    // run active model components
    // first diffuse the hillslopes. Currently two options. Perhaps a flag could be
    // added so that we don't have to keep adding if statements as more
    // hillslope rules are added?
    if (hillslope)
    {
      if (nonlinear)
      {
        //soil_diffusion_fv_nonlinear();
        MuddPILE_nl_soil_diffusion_nouplift();
      }
      else
      {
        soil_diffusion_fd_linear();
      }
    }

    // sediment will have moved into channels. This is assumed to
    // be immediately removed by fluvial processes, so we use the
    // 'wash out' member function
    wash_out();

    // now for fluvial erosion
    if (fluvial)
    {
      // currently the only option is stream power using the FASTSCAPE
      // algorithm so there are no choices here
      fluvial_incision_with_uplift();
    }

    // now for isostacy. As of 17/6/2014 this wasn't working
    if (isostasy)
    {
      if (flexure)
        //flexural_isostasy( 0.00001 );
        flexural_isostasy_alt( );
      else
        Airy_isostasy();
    }

    // this writes a file if it is on the print interval, and writes some erosion
    // data to an array always so that steady state may be tested
    write_report();

    //update the time
    current_time += timeStep;

    // write at every print interval
    if ( print_interval > 0 && (print % print_interval) == 0)
    {
      print_rasters_and_csv( frame );
      ++frame;
    }
    if (quiet == false) cout << "\rTime: " << current_time << " years" << flush;
    ++print;

    // check to see if steady state has been achieved
    check_steady_state();

  } while (check_end_condition() == false);

  if ( print_interval == 0 || (print_interval > 0 && ((print-1) % print_interval) != 0))
  {
    // if the print interval is 0, we set the frame to zero
    if (print_interval == 0)
    {
      frame = 0;
    }
    print_rasters_and_csv( frame );
  }

  // reset the current frame
  current_frame = frame;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is a wrapper function used to drive a model run
// it checks parameters and flags from the data members
//
// Similar to run_componets but instead of running fluvial
// and uplift in seperate steps it inserts the fluvial
// and uplift matrices into the nonlinear hillslope
// solver.
//
// This version allows variable K and U rasters to be used.
// You can also switch on an adaptive timestep
//
// SMM, 3/09/2017
//  updated 6/9/2017 with adaptive timestep
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::run_components_combined( LSDRaster& URaster, LSDRaster& KRaster, bool use_adaptive_timestep)
{
  //recording = false;
  cycle_number = 1;
  total_erosion = 0;
  max_erosion = 0;
  min_erosion = -99;
  switch_delay = 0;
  time_delay = 0;

  stringstream ss, ss_root;

  // some fields for uplift and fluvial incision
  Array2D<float> uplift_field;
  Array2D<float> fluvial_incision_rate_field;

  // set the frame to the current frame plus 1
  int frame = current_frame;
  do
  {
    // Record current topography
    zeta_old = RasterData.copy();

    // run active model components
    // first diffuse the hillslopes. Currently two options. Perhaps a flag could be
    // added so that we don't have to keep adding if statements as more
    // hillslope rules are added?
    if (hillslope)
    {
      if (nonlinear)
      {
        //soil_diffusion_fv_nonlinear();
        MuddPILE_nl_soil_diffusion_nouplift();
      }
      else
      {
        soil_diffusion_fd_linear();
      }
    }

    // sediment will have moved into channels. This is assumed to
    // be immediately removed by fluvial processes, so we use the
    // 'wash out' member function
    wash_out();

    // now for fluvial erosion
    if (fluvial)
    {
      // currently the only option is stream power using the FASTSCAPE
      // algorithm so there are no choices here

      if (use_adaptive_timestep)
      {
        //cout << "I am using an adaptive timestep" << endl;
        fluvial_incision_with_variable_uplift_and_variable_K_adaptive_timestep( URaster, KRaster );
      }
      else
      {
        fluvial_incision_with_variable_uplift_and_variable_K( URaster, KRaster );
      }
    }

    //update the time
    current_time += timeStep;

    // now see if the time has exceeded the next print time
    if (current_time > next_printing_time)
    {
      do
      {
        next_printing_time+=float_print_interval;
      } while(next_printing_time < current_time);
      print_rasters_and_csv( frame );
      ++frame;
    }
    if (quiet == false) cout << "\rTime: " << current_time << " years" << flush;

  } while (check_end_condition() == false);

  // reset the current frame
  current_frame = frame;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is a wrapper function used to drive a model run
// it checks parameters and flags from the data members
//
// This includes a number of columns for tracking erosion rates
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::run_components_combined_cell_tracker( vector<LSDParticleColumn>& CRNColumns,
                      vector<LSDParticleColumn>& eroded_cells,
                      int startType, double startDepth, double particle_spacing,
                      LSDCRNParameters& CRNParams)
{

  // first you need to get the locations of the columns
  int N_pcolumns = CRNColumns.size();
  vector<double> row_vec;
  vector<double> col_vec;

  vector<LSDParticleColumn> e_cells(N_pcolumns);

  for (int i = 0; i<N_pcolumns; i++)
  {
    row_vec.push_back( CRNColumns[i].getRow() );
    col_vec.push_back( CRNColumns[i].getCol() );
  }

  stringstream ss, ss_root;


  //cout << "\n\nLine 2604, data[10][10]: " << RasterData[10][10] << endl;

  // set the frame to the current frame
  int frame = current_frame;
  int print = 1;
  do
  {
    // Record current topography
    zeta_old = RasterData.copy();

    // run active model components
    // first diffuse the hillslopes. Currently two options. Perhaps a flag could be
    // added so that we don't have to keep adding if statements as more
    // hillslope rules are added?
    if (hillslope)
    {
      //cout << "Time: " << current_time << " nonlinear: " << nonlinear << endl;
      if (nonlinear)
      {
        //soil_diffusion_fv_nonlinear();
        MuddPILE_nl_soil_diffusion_nouplift();
      }
      else
        soil_diffusion_fd_linear();
    }

    //cout << "Line 2600, data[10][10]: " << RasterData[10][10] << endl;

    // sediment will have moved into channels. This is assumed to
    // be immediately removed by fluvial processes, so we use the
    // 'wash out' member function
    wash_out();

    // now for fluvial erosion
    if (fluvial)
    {
      // currently the only option is stream power using the FASTSCAPE
      // algorithm so there are no choices here
      fluvial_incision_with_uplift();
    }

    //update the time
    current_time += timeStep;

    //cout << "/n/nLine 2604, data[10][10]: " << RasterData[10][10] << endl;

    // get the erosion rates at the cells
    for (int i = 0; i<N_pcolumns; i++)
    {
      int row, col;
      row = row_vec[i];
      col = col_vec[i];

      double startxLoc = double(CRNColumns[i].getCol())*DataResolution+XMinimum;
      double startyLoc = double(NRows-CRNColumns[i].getRow()-1)*DataResolution+YMinimum;

      double this_zeta_old = zeta_old[row][col];
      double this_zeta_new = RasterData[row][col];



      //double test = RasterData[row][col];
      //cout << "zo " << this_zeta_old << " RD: " <<  RasterData[row][col] << endl;
      //test = 5;
      //cout << "test is: " << test << " RD: " <<  RasterData[row][col] << endl;


      double this_uplift_rate = get_uplift_rate_at_cell(row,col);

      //cout << "\n\n\nCRN at [" << row << "]["<<col<<"], znew: " << this_zeta_new
      //     << " zold: " << this_zeta_old << " uplift rate: " << this_uplift_rate << endl;
      //CRNColumns[i].print_particle_properties_to_screen(CRNParams);


      LSDParticleColumn this_eroded_column =
                 CRNColumns[i].update_CRN_list_rock_only_eros_limit_3CRN(
                       timeStep, this_uplift_rate,
                       startType, startDepth,
                       startxLoc,  startyLoc,
                       this_zeta_old,this_zeta_new,
                       particle_spacing, CRNParams);
      e_cells[i] = this_eroded_column;

      //cout << "\n\n\nAFTER CHANGE!!! CRN at [" << row << "]["<<col<<"], znew: " << this_zeta_new
      //     << " zold: " << this_zeta_old << " uplift rate: " << this_uplift_rate << endl;
      //CRNColumns[i].print_particle_properties_to_screen(CRNParams);

    }

    //CRNColumns[1].print_particle_properties_to_screen(CRNParams);

    //cout << "Line 2643, data[10][10]: " << RasterData[10][10] << endl;

    // write at every print interval
    if ( print_interval > 0 && (print % print_interval) == 0)
    {
      // calculate erosion rate and place it in the erosion data member for
      // raster printing
      erosion = get_total_erosion_rate_over_timestep();
      print_rasters( frame );

      print_average_erosion_and_apparent_erosion( frame, CRNColumns, CRNParams);
      //print_column_erosion_and_apparent_erosion( frame, CRNColumns, CRNParams);

      ++frame;
    }
    if ( quiet == false) cout << "\rTime: " << current_time << " years" << flush;
    ++print;

    //cout << "Line 2693, data[10][10]: " << RasterData[10][10] << endl;

    // check to see if steady state has been achieved
    //check_steady_state();

  }  while ( check_end_condition() == false);

  eroded_cells = e_cells;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function intiates a vector of LSDParticleColumns at steady state
// the zeta values are distrubuted over the surface
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<LSDParticleColumn> LSDRasterModel::initiate_steady_CRN_columns(int column_spacing,
                    vector<int>& CRNcol_rows, vector<int>& CRNcol_cols,
                    double rho_r, double this_U, int startType, double startDepth,
                    double particle_spacing, LSDCRNParameters& CRNParam)
{
  cout << "Initialising some CRN columns." << endl;

  vector<LSDParticleColumn> CRNParticleColumns_vec;

  // get the number of particle rows and columns
  // note there are less rows since these are boundary nodes that
  // have zero erosion and are not included
  int NPRows = (NRows-2)/column_spacing+1;
  int NPCols = NCols/column_spacing+1;
  cout << "Rows: " << NRows << " and NPRows: " << NPRows << endl;
  cout << "Cols: " << NCols << " and NPCols: " << NPCols << endl;
  vector<int> rows_vec;
  vector<int> cols_vec;

  int this_row, this_col;
  double this_x,this_y;
  double zeta_at_cell;
  double eff_U = rho_r*this_U/10;     // this converts the uplift to effective
                                      // uplift in units g/cm^2/yr,
                                      // required by the cosmo particles

  // loop through the columns, creating steady particles
  for (int prow = 0; prow<NPRows;  prow++)
  {
    this_row = prow*column_spacing+1;
    for (int pcol = 0; pcol<NPCols; pcol++)
    {
      this_col = pcol*column_spacing;
      // make sure you dont go out of bounds
      if(this_row >= NRows-1)
      {
        cout << "Danger, you were about to go out of the raster domain when making a CRN column\n";
        this_row = NRows-2;
      }
      if(this_col >= NCols)
      {
        cout << "Danger, you were about to go out of the raster domain when making a CRN column\n";
        this_col = NCols-1;
      }

      // push back the row and col vec
      rows_vec.push_back(this_row);
      cols_vec.push_back(this_col);

      // get the x and y locations
      this_x = double(this_col)*DataResolution+XMinimum;
      this_y = double(NRows-this_row-1)*DataResolution+YMinimum;

      // initiate an empty column
      LSDParticleColumn temp_col;
      temp_col.set_Row_and_Col(this_row,this_col);

      // get the elevation of the surface
      zeta_at_cell = float(RasterData[this_row][this_col]);

      // make a steady state column
      //cout << "LINE 3032, initiating SS column" << endl;
      temp_col.initiate_SS_cosmo_column_3CRN(startType, this_x, this_y, startDepth, particle_spacing,
                                             zeta_at_cell, eff_U, CRNParam);
      //cout << "LINE 3035, Initiated column" << endl;

      // make sure this is a rock only simulation
      double ST = 0;
      temp_col.set_SoilThickness(ST);
      temp_col.set_RockDensity(rho_r);

      // add the column to the vector
      CRNParticleColumns_vec.push_back(temp_col);
      //cout << "Got col, row: " << this_row << " and col: " << this_col << endl;
    }
  }

  //cout << "/n/nLine 2788, data[10][10]: " << RasterData[10][10] << endl;
  CRNcol_rows = rows_vec;
  CRNcol_cols = cols_vec;
  return CRNParticleColumns_vec;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function simply calls the run_components function
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::run_model( void )
{
  // file to write details of each time step
  // This will be opened by the write_report method
  int run = 1;

  // loop through the number of runs. The number of runs is stored as a data member
  do
  {
    //initial_steady_state = false;      // SMM: not clear why this is set to false
    // SMM: I have turned this off because it was overriding the periodic forcing

    recording = false;                 // SMM: not clear why this is set to false

    if ( initialized == false &&  quiet == false)
    {
      cout << "Model has not been initialized with a parameter file." << endl;
      cout << "All values used are defaults" << endl;
    }

    current_time = 0;
    run_components();
    ++run;
  } while (run <= num_runs);

  final_report();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function starts from a steady state.
// apparently the steady state data is stored as a data member.
// SMM: not sure if this is the best thing to do since it adds a lot
// of memory overhead. Could this be done by using a steady state raster
// as an argument?
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::run_model_from_steady_state( void )
{
  // file to write details of each time step
  // This will be opened by the write_report method
  RasterData = steady_state_data;
  reset_model();

  if ( initialized == false &&  quiet == false)
  {
    cout << "Model has not been initialized with a parameter file." << endl;
    cout << "All values used are defaults" << endl;
  }
  if ( initial_steady_state == false)
  {
    cout << "Model has not been set to steady state yet" << endl;
    cout << "Run LSDRasterModel::reach_steady_state( float tolerance ) first" << endl;
  }
  // Generate random noise

  int run = 1;
  stringstream ss, ss_root;
  do{
    current_time = 0;
    run_components();
    ++run;
  } while (run <= num_runs);

  // If the last output wasn't written, write it

  if ( quiet == false) cout << "\nModel finished!\n" << endl;
  final_report();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function brings a model run to steady state
// It does not simply run at a constant erosion rate, but rather
// runs using periodic forcing. Tests have shown this can greatly reduce
// the time required to reach steady state
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::reach_steady_state( void )
{
  // file to write details of each time step
  // This will be opened by the write_report method
  initial_steady_state = false;        // this means that the model
                                       // knows SS has not been reached so it
                                       // will not 'record' model results
  current_time = 0;
  total_erosion = 0;
  max_erosion = 0;
  min_erosion = -99;

  // these 'swap' functions are here because the run to steady state
  // uses default paramaters but the original model parameters are
  // stored in these 'swap' data members and restored after running
  // the model to steady state.
  int K_mode_swap = K_mode;
  int D_mode_swap = D_mode;
  float K_amp_swap = K_amplitude;
  float endTime_swap = endTime;
  float period_swap = periodicity;
  int period_mode_swap = period_mode;
  float print_interval_swap = print_interval;
  bool reporting_swap = reporting;
  string name_swap = name;

  if ( initialized == false &&  quiet == false)
  {
    cout << "Model has not been initialized with a parameter file." << endl;
    cout << "All values used are defaults" << endl;
  }

  // Generate random noise
  random_surface_noise(0, noise);
  // Fill the topography
  LSDRaster *temp;
  temp = new LSDRaster(*this);
  float thresh_slope = 0.00001;
  *temp = fill(thresh_slope);
  RasterData = temp->get_RasterData();
  delete temp;

  // Run model with some modest fluvial forcing
  K_mode = 1;                    // THis means a sinusoidal variation in
                                 // channel forcing
  K_amplitude = K_fluv * 0.3;
  endTime = 0;                   // SMM: not sure why this is set to zero.
  cout << "Running to steady state, periodicity is: " << periodicity << endl;
  //periodicity = 2000;
  period_mode = 1;
  cycle_steady_check = true;       // this means that the steady state check
                                   // is run based on erosion from one cycle
                                   // to the next
  print_interval = 0;
  reporting = false;

  if ( quiet == false)
  {
    cout << "Producing steady state profile" << endl;
  }
  // because the cycle_steady_check is true, it means that this will run
  // through several periods until the surface does not change between cycles

  // switch the end time mode to be periodic
  short endTimeModeSwap = endTime_mode;
  endTime_mode = 2;

  // run components
  run_components();

  if ( quiet == false)
  {
    cout << "Finished producing cyclic steady" << endl;
  }

  // switch back to the original endTime_mode
  endTime_mode =  endTimeModeSwap;

  // Now run with static forcing
  K_mode = K_mode_swap;
  D_mode = D_mode_swap;
  K_amplitude = K_amp_swap;
  endTime = timeStep*10;
  cycle_steady_check = false;
  initial_steady_state = false;
  current_time = 0;

  if ( quiet == false) cout << "Did cyclic steady, now running at constant forcing for a while." << endl;
  cout << "Timestep is: " << timeStep << endl;
  run_components();
  if ( quiet == false) cout << "Forced from constant steady state, exiting steady state routine." << endl;

  endTime = endTime_swap;
  periodicity = period_swap;
  period_mode = period_mode_swap;
  cycle_steady_check = false;
  print_interval = print_interval_swap;
  reporting = reporting_swap;

  steady_state_data = Array2D<float> (NRows, NCols, 0.0);
  steady_state_data = RasterData;

  // set the initial steady state flag to true
  initial_steady_state = true;
  cout << "Line 2526 Finished run_to_steady_state, initial_steady_state: " << initial_steady_state << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


void LSDRasterModel::soil_diffusion_fv( void )
{
  static bool defined = false;

  static int problem_dimension;
  static float inv_dx_S_c_squared, inv_dy_S_c_squared, dx_front_term, dy_front_term;
  static vector<int> vec_k_value_i_j;
  static vector<int> vec_k_value_ip1_j;
  static vector<int> vec_k_value_im1_j;
  static vector<int> vec_k_value_i_jp1;
  static vector<int> vec_k_value_i_jm1;
  static float iteration_tolerance = 0.01;

  Array2D <float> fluvial_temp(NRows, NCols, 0.0);
  float south, north;

  if (boundary_conditions[2][0] == 'b')
  {
    south = 0.0;
    north = current_time*get_max_uplift();

    cout << south << ", " << north << endl;
  }
  else
  {
    cout << "Model currently not built to cope with hillslope diffusion using these boundary conditions" << endl;
    cout << "Feature implementation required" << endl;
    exit(1);
  }


  if (defined == false)
  {
  mtl_initiate_assembler_matrix(problem_dimension, inv_dx_S_c_squared, inv_dy_S_c_squared,
                  dx_front_term, dy_front_term, vec_k_value_i_j, vec_k_value_ip1_j,
      vec_k_value_im1_j, vec_k_value_i_jp1, vec_k_value_i_jm1);
  }
  defined = true;

  nonlinear_creep_timestep(fluvial_temp, iteration_tolerance, problem_dimension,
      inv_dx_S_c_squared, inv_dy_S_c_squared,  dx_front_term, dy_front_term, vec_k_value_i_j,
      vec_k_value_ip1_j, vec_k_value_im1_j, vec_k_value_i_jp1, vec_k_value_i_jm1,
      south, north);
}




mtl::compressed2D<float> LSDRasterModel::generate_fd_matrix( int dimension, int size, bool periodic )
{
  int num_neighbours, num_neighbours_;
  int row, col;
  float r = get_D() * timeStep / (DataResolution * DataResolution);
  float r_ = get_D() * timeStep / pow(DataResolution * 1.4142135623, 2);
  int width, height;

  if (dimension == 0)    // North - south
  {
    width = NCols;
    height = NRows - 2;
  }
  else          // East - west
  {
    width = NCols - 2;
    height = NRows;
  }

  mtl::compressed2D<float> matrix(size, size);
  matrix = 0.0;
  mtl::mat::inserter< mtl::compressed2D<float> > ins(matrix);

  for (int i=0; i<size; ++i)
  {
    row = i / width;
    col = i % width;
    num_neighbours = 4;
    num_neighbours_ = 4;
    // left
    if (col > 0)
    {
      ins[i][i-1] << -r;
    }
    else if (dimension == 0)
    {
      if ( periodic == false)
        --num_neighbours;
      else
        ins[i][i+width-1] << -r;
    }

    // right
    if (col < width - 1)
    {
      ins[i][i+1] << -r;
    }
    else if (dimension == 0)
    {
      if ( periodic == false)
        --num_neighbours;
      else
        ins[i][i-width+1] << -r;
    }

    // up
    if (row > 0)
    {
      ins[i][i-width] << -r;
    }
    else if (dimension == 1)
    {
      if ( periodic == false)
        --num_neighbours;
      else
        ins[i][i+(width*(NCols-1))] << -r;
    }

    // down
    if (row < height-1)
    {
      ins[i][i+width] << -r;
    }
    else if (dimension == 1)
    {
      if (periodic == false)
        --num_neighbours;
      else
        ins[i][i-(width*(NCols-1))] << -r;
    }

    // Diagonals
    // Upper left
    if (row > 0 && col > 0)
      ins[i][i-width-1] << -r_;
    else if (dimension == 0 && row > 0)
      if (periodic == false)
        --num_neighbours_;
      else{
        ins[i][i-1] << -r_;}
    else if (dimension == 1 && col > 0)
    {
      if (periodic == false)
        --num_neighbours_;
      else
        ins[i][i+(width*(NCols-1))-1] << -r;
    }


    // Upper right
    if (row > 0 && col < width-1)
      ins[i][i-width+1] << -r_;
    else if (dimension == 0 && row > 0)
    {
      if (periodic == false)
        --num_neighbours_;
      else
        ins[i][i-(2*width)+1] << -r_;
    }
    else if (dimension == 1 && col < width-1)
    {
      if (periodic == false)
        --num_neighbours_;
      else
        ins[i][i+(width*(NCols-1))+1] << -r_;
    }

    // Lower left
    if (row < height-1 && col > 0)
      ins[i][i+width-1] << -r_;
    else if (dimension == 0 && row < height-1)
    {
      if (periodic == false)
        --num_neighbours_;
      else
        ins[i][i+(2*width)-1] << -r_;
    }

    else if (dimension == 1 && col > 0)
    {
      if (periodic == false)
        --num_neighbours_;
      else
        ins[i][col-1] << -r_;
    }

    // Lower right
    if (row < height-1 && col < width-1)
    {
      ins[i][i+width+1] << -r_;
    }
    else if (dimension == 0 && row < height-1)
    {
      if (periodic == false)
        --num_neighbours_;
      else
        ins[i][i+1] << -r_;
    }

    else if (dimension == 1 && col < width-1)
    {
      if (periodic == false)
        --num_neighbours_;
      else
        ins[i][col+1] << -r_;
    }

    ins[i][i] << num_neighbours*r + 1 + num_neighbours_ * r_;
  }
  return matrix;
}

mtl::dense_vector <float> LSDRasterModel::build_fd_vector(int dimension, int size)
{
  int vector_pos = 0;
  mtl::dense_vector <float> data_vector(size);
  float push_val;
  float r = get_D() * timeStep / (DataResolution * DataResolution);
  int start_i, end_i;
  int start_j, end_j;

  if (dimension == 0)
  {
    start_i = 1; end_i = NRows-2;
    start_j = 0; end_j = NCols-1;
  }
  else       // East - west
  {
    start_i = 0; end_i = NRows-1;
    start_j = 1; end_j = NCols-2;
  }

  for (int i=start_i; i<=end_i; ++i)
  {
    for (int j=start_j; j<=end_j; ++j)
    {
      push_val = RasterData[i][j];
      if (dimension == 0)
      {
        if (i==1)
          push_val += RasterData[0][j] * r;
        else if (j==NRows-2)
          push_val += RasterData[NRows-1][j] * r;
      }

      else if (dimension == 1)
      {
      if (j == 1)
        push_val += RasterData[i][0] * r;
      else if (j == NCols-2)
        push_val += RasterData[i][NCols-1] * r;
      }

      data_vector[vector_pos] = push_val;
      ++vector_pos;
    }
  }

  return data_vector;
}

/*
void LSDRasterModel::repack_fd_vector(mtl::dense_vector <float> &data_vector, int dimension)
{
  int vector_pos = 0;
  if (dimension == 1)       // East - west
  {
    for (int i=0; i<NRows; ++i)
    {
      for (int j=1; j<NCols-1; ++j)
      {
        RasterData[i][j] = data_vector[vector_pos];
        ++vector_pos;
      }
    }
  }
  if (dimension == 0)       // North - south
  {
    for (int j=0; j<NCols; ++j)
    {
      for (int i=1; i<NRows-1; ++i)
      {
        RasterData[i][j] = data_vector[vector_pos];
        ++vector_pos;
      }
    }
  }
}
*/

void LSDRasterModel::soil_diffusion_fd_linear( void )
{
  short dimension;
  bool periodic;
  int size;

  interpret_boundary(dimension, periodic, size);
  //cout << "Periodic " << periodic << endl;

  mtl::compressed2D <float> matrix = generate_fd_matrix(dimension, size, periodic);
  // Unpack data
  mtl::dense_vector <float> data_vector = build_fd_vector(dimension, size);

  if ( quiet == false && name == "debug" && size < 100)
  {
    cout << "Data: " << endl;
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        cout << RasterData[i][j] << " ";
      }
      cout << endl;
    }
    cout << "Matrix: " << endl;
    for (int i=0; i<size; ++i)
    {
      for (int j=0; j<size; ++j)
      {
        cout << matrix[i][j] << " ";
      }
      cout << endl;
    }
    cout << "Vector: " << endl;
    for (int i=0; i<size; ++i)
      cout << data_vector[i] << endl;
  }

  mtl::dense_vector <float> output(size);
  // Set up preconditioner
  itl::pc::ilu_0 < mtl::compressed2D<float> > P(matrix);
  // Iterator condition
  itl::basic_iteration <float> iter(data_vector, 200, 1e-6);
  // Matrix solver
  itl::bicgstab(matrix, output, data_vector, P, iter);


  repack_vector(output, dimension);
  if (quiet == false && name == "debug" && size <100)
  {
    cout << "Output: " << endl;
    for (int i=0; i<size; ++i)
      cout << output[i] << endl;
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        cout << RasterData[i][j] << " ";
      }
      cout << endl;
    }
  }
}

mtl::compressed2D<float> LSDRasterModel::generate_fv_matrix( int dimension, int size, bool periodic )
{
  float A, B, C, D;
  float front = timeStep * get_D() / (DataResolution*DataResolution);
  float inv_term = 1 / (DataResolution * DataResolution * S_c * S_c);

  int p = 0;  // positioner for matrix insertion
  int offset;
  int start_i, start_j, end_i, end_j;
  mtl::compressed2D <float> matrix(size, size);
  matrix = 0.0;
  mtl::mat::inserter< mtl::compressed2D<float> > ins(matrix);

  if (dimension == 0)
  {
    start_i = 1; end_i = NRows-2;
    start_j = 0; end_j = NCols-1;
    offset = NCols;
  }
  else
  {
    start_i = 0; end_i = NRows-1;
    start_j = 1; end_j = NCols-2;
    offset = NCols - 2;
  }

  for (int i=start_i; i<=end_i; ++i)
  {
    for (int j=start_j; j<=end_j; ++j)
    {
      A = (i==0)    ? 0 : front / (1 - pow(RasterData[i][j] - RasterData[i-1][j], 2) * inv_term);
      B = (j==NCols-1)  ? 0 : front / (1 - pow(RasterData[i][j] - RasterData[i][j+1], 2) * inv_term);
      C = (i==NRows-1)   ? 0 : front / (1 - pow(RasterData[i][j] - RasterData[i+1][j], 2) * inv_term);
      D = (j==0)    ? 0 : front / (1 - pow(RasterData[i][j] - RasterData[i][j-1], 2) * inv_term);

      if (periodic)
      {
        if (i==0)     A = front / (1 - pow(RasterData[i][j] - RasterData[NRows-1][j], 2) * inv_term);
        else if (j==NCols-1)   B = front / (1 - pow(RasterData[i][j] - RasterData[i][0], 2) * inv_term);
        else if (i==NRows-1)   C = front / (1 - pow(RasterData[i][j] - RasterData[0][j], 2) * inv_term);
        else if (j==0)     D = front / (1 - pow(RasterData[i][j] - RasterData[i][NCols-1], 2) * inv_term);

      }

      ins[p][p] << 1 + A + B + C + D;
      if (j != start_j)
        ins[p][p-1] << -D;
      else if (periodic && dimension == 0 )
        ins[p][p+offset-1] << -D;
      if (j != end_j)
        ins[p][p+1] << -B;
      else if (periodic && dimension == 0 )
        ins[p][p-offset+1] << -B;
      if (i != start_i)
        ins[p][p-offset] << -A;
      else if (periodic && dimension == 1 )
        ins[p][p+(offset*(NCols-1))] << -A;
      if (i != end_i)
        ins[p][p+offset] << -C;
      else if (periodic && dimension == 1 )
        ins[p][p-(offset*(NCols-1))] << -C;

      ++p;
    }
  }
  return matrix;
}

mtl::dense_vector <float> LSDRasterModel::build_fv_vector( int dimension, int size )
{
  float front = timeStep * get_D() / (DataResolution*DataResolution);
  float inv_term = 1 / (DataResolution * DataResolution * S_c * S_c);

  mtl::dense_vector <float> data_vector(size);
  int p = 0;    // vector positioner
  int start_i, end_i;
  int start_j, end_j;
  float push_val;

  if (dimension == 0)
  {
    start_i = 1; end_i = NRows-2;
    start_j = 0; end_j = NCols-1;
  }
  else
  {
    start_i = 0; end_i = NRows-1;
    start_j = 1; end_j = NCols-2;
  }

  for (int i=start_i; i<=end_i; ++i)
  {
    for (int j=start_j; j<=end_j; ++j)
    {
      push_val = zeta_old[i][j];
//      cout << push_val << endl;

      if (dimension == 0)
      {
        if (i==1)
          push_val += zeta_old[0][j] * front /
            (1 - pow(RasterData[i][j]-RasterData[0][j],2) * inv_term);
        if (i==NRows-2)
          push_val += zeta_old[NRows-1][j] * front /
            (1 - pow(RasterData[i][j]-RasterData[NRows-1][j],2) * inv_term);
      }
      else if (dimension == 1)
      {
        if (j==1)
          push_val += zeta_old[i][0] * front /
            (1 - pow(RasterData[i][j]-RasterData[i][0],2) * inv_term);
        if (j==NCols-2)
          push_val += zeta_old[i][NCols-1] * front /
            (1 - pow(RasterData[i][j]-RasterData[i][NCols-1],2) * inv_term);
      }

      data_vector[p] = push_val;
      ++p;
    }
  }
  return data_vector;
}

void LSDRasterModel::repack_vector(mtl::dense_vector <float> &data_vector, int dimension)
{
  int start_i, end_i;
  int start_j, end_j;
  int p = 0;

  if (dimension == 0)
  {
    start_i = 1; end_i = NRows-2;
    start_j = 0; end_j = NCols-1;
  }
  else
  {
    start_i = 0; end_i = NRows-1;
    start_j = 1; end_j = NCols-2;
  }

  for (int i = start_i; i<=end_i; ++i)
  {
    for (int j = start_j; j<=end_j; ++j)
    {
      RasterData[i][j] = data_vector[p];
      ++p;
    }
  }
}


void LSDRasterModel::soil_diffusion_fv_nonlinear( void )
{
  int max_iter = 200, iter = 0;
  float epsilon = 0.00001;
  float max_diff;
  Array2D <float> last_iteration;

  short dimension = 0;
  bool periodic;
  int size;
  interpret_boundary(dimension, periodic, size);

  do {
  last_iteration = RasterData.copy();
  // A
  mtl::compressed2D <float> matrix = generate_fv_matrix(dimension, size, periodic);
  // b
  mtl::dense_vector <float> data_vector = build_fv_vector(dimension, size);

  if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
  {
  cout << "Data: " << endl;
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      cout << RasterData[i][j] << " ";
    }
    cout << endl;
  }
  cout << "Matrix: " << endl;
  for (int i=0; i<size; ++i)
  {
    for (int j = 0; j<size; ++j)
    {
      cout <<matrix[i][j] << " ";
    }
    cout << endl;
  }
  cout << "Vector: " << endl;
  for (int i=0; i<size; ++i)
    cout << data_vector[i] << endl;
  }

  // x
  mtl::dense_vector <float> output(size);
  // Set up preconditioner
  itl::pc::ilu_0 < mtl::compressed2D<float> > P(matrix);
  // Iterator condition
  itl::basic_iteration <float> iter(data_vector, 200, 1e-6);
  // Matrix solver
  itl::bicgstab(matrix, output, data_vector, P, iter);

  repack_vector(output, dimension);
  /*
  if (NRows <= 10 && NCols <= 10)
  {
  cout << "Output: " << endl;
  for (int i=0; i<size; ++i)
    cout << output[i] << endl;
  cout << "New data: " << endl;
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      cout << RasterData[i][j] << " ";
    }
    cout << endl;
  }
  }
  */

    max_diff = 0;
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        if (abs(RasterData[i][j] - last_iteration[i][j]) > max_diff)
          max_diff = RasterData[i][j] - last_iteration[i][j];
      }
    }
  if (quiet == false && name == "debug" && size <100)
  {
    cout << "Output: " << endl;
    for (int i=0; i<size; ++i)
      cout << output[i] << endl;
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        cout << RasterData[i][j] << " ";
      }
      cout << endl;
    }
  }
  ++iter;
  } while (max_diff > epsilon && iter < max_iter);

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This takes the model and calculates the steady state fluvial surface derived from
// a chi map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::fluvial_snap_to_steady_state(float U)
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // need to fill the raster to ensure there are no internal base level nodes
  float slope_for_fill = 0.0001;
  LSDRaster filled = temp.fill(slope_for_fill);
  LSDFlowInfo flow(boundary_conditions, filled);

  float K = get_K();
  float m_exp = get_m();
  float n_exp = get_n();
  float area_threshold = 0;
  float m_over_n = m_exp/n_exp;
  float A_0 = 1;
  float one_over_n = 1/n_exp;

  LSDRaster ChiValues = flow.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0,area_threshold);
  float thisChi;
  // now we calculate the elevations assuming that the elevation at chi = 0 is 0
  // This is based on equation 4 from mudd et al 2014 JGR-ES
  for (int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      thisChi = ChiValues.get_data_element(row,col);
      if (thisChi == NoDataValue)
      {
        zeta[row][col] = NoDataValue;
      }
      else
      {
        if (n_exp == 1)
        {
          zeta[row][col] = (U/K)*thisChi;
        }
        else
        {
          zeta[row][col] = pow( (U/K),one_over_n)*thisChi;
        }

      }
    }
  }

  RasterData = zeta.copy();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This snaps to steady based on an input file with elevations and node indicies
// overloaded from the previous function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::fluvial_snap_to_steady_variable_K_variable_U(LSDRaster& K_values, 
                                  LSDRaster& U_values, LSDSpatialCSVReader& source_points_data, 
                                  bool carve_before_fill)
{

  string column_name = "elevation(m)";


  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // need to fill the raster to ensure there are no internal base level nodes
  cout << "I am going to carve and fill" << endl;
  float slope_for_fill = 0.0001; 
  LSDRaster filled_topography,carved_topography;
  if(carve_before_fill)
  {
    cout << "Carving and filling." << endl;
    carved_topography = temp.Breaching_Lindsay2016();
    filled_topography = carved_topography.fill(slope_for_fill);
  }
  else
  {
    cout << "Filling." << endl;
    filled_topography = temp.fill(slope_for_fill);
  }
  cout << "Getting the flow info. This might take some time." << endl;
  LSDFlowInfo flow(boundary_conditions, filled_topography);
  // update the raster
  zeta = filled_topography.get_RasterData();

  // Get the local node index as well as the elevations
  vector<int> ni = source_points_data.get_nodeindices_from_lat_long(flow);
  vector<float> elev = source_points_data.data_column_to_float(column_name);
  // make the map
  cout << "I am making an elevation map. This will not work if points in the raster lie outside of the csv channel points." << endl;
  map<int,float> elevation_map;
  for(int i = 0; i< int(ni.size()); i++)
  {
    elevation_map[ ni[i] ] =  elev[i]; 
  }

  float m_exp = get_m();
  float n_exp = get_n();
  float one_over_n = 1/n_exp;

  //float FP_NDV = K_values.get_NoDataValue();
  float FP_value, K_value, U_value, receiver_elev, parenth_term, area_pow;


  // Step one, create donor "stack" etc. via FlowInfo
  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx;

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;

  // Step two calculate new height
  //for (int i=numNodes-1; i>=0; --i)
  for (int i=0; i<numNodes; ++i)
  {

    // get the information about node relationships from the flow info object
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DR2;


    // check if this is a baselevel node
    if(node == receiver)
    {
      //cout << "This is a base level node. I don't update this node." << endl;
    }
    else if(elevation_map.find(node) != elevation_map.end())
    {
      //cout << "This is one of the fixed channel nodes. I don't update this node." << endl;
      FP_value = elevation_map[node];
      zeta[row][col]= FP_value;
    }
    else
    {
      // Get the data from the individual points
      K_value = K_values.get_data_element(row,col);
      U_value = U_values.get_data_element(row,col);

      // get the distance between nodes. Depends on flow direction
      switch (flow.retrieve_flow_length_code_of_node(node))
      {
        case 0:
          dx = -99;
        break;
        case 1:
          dx = DataResolution;
        break;
        case 2:
          dx = dx_root2;
        break;
        default:
          dx = -99;
        break;
      }


      receiver_elev = zeta[receiver_row][receiver_col];
      area_pow = pow(drainageArea,m_exp);
      parenth_term = U_value/(K_value*area_pow);
      zeta[row][col] = receiver_elev+ dx*( pow(parenth_term,one_over_n));

    }
    



  }
  //return LSDRasterModel(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This snaps to steady based on an input file with elevations and node indicies
// overloaded from the previous function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::fluvial_snap_to_steady_variable_K_variable_U(LSDRaster& K_values, LSDRaster& U_values, string csv_of_fixed_channel, bool carve_before_fill)
{
  LSDRasterInfo RI(K_values);

  cout << "I am reading points from the file: "+ csv_of_fixed_channel << endl;
  LSDSpatialCSVReader source_points_data( RI,csv_of_fixed_channel );

  //*****************************
  // Working from here
  // Need to get out map that has key as node index and value as elevation
  //******************************
  string column_name = "elevation(m)";


  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // need to fill the raster to ensure there are no internal base level nodes
  float slope_for_fill = 0.0001; 
  LSDRaster filled_topography,carved_topography;
  if(carve_before_fill)
  {
    cout << "Carving and filling." << endl;
    carved_topography = temp.Breaching_Lindsay2016();
    filled_topography = carved_topography.fill(slope_for_fill);
  }
  else
  {
    cout << "Filling." << endl;
    filled_topography = temp.fill(slope_for_fill);
  }
  cout << "Getting the flow info. This might take some time." << endl;
  LSDFlowInfo flow(boundary_conditions, filled_topography);
  // update the raster
  zeta = filled_topography.get_RasterData();


  // Get the local node index as well as the elevations
  vector<int> ni = source_points_data.get_nodeindices_from_lat_long(flow);
  vector<float> elev = source_points_data.data_column_to_float(column_name);
  // make the map
  cout << "I am making an elevation map. This will not work if points in the raster lie outside of the csv channel points." << endl;
  map<int,float> elevation_map;
  for(int i = 0; i< int(ni.size()); i++)
  {
    elevation_map[ ni[i] ] =  elev[i]; 
  }

  float m_exp = get_m();
  float n_exp = get_n();
  float one_over_n = 1/n_exp;

  //float FP_NDV = K_values.get_NoDataValue();
  float FP_value, K_value, U_value, receiver_elev, parenth_term, area_pow;


  // Step one, create donor "stack" etc. via FlowInfo
  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx;

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;

  // Step two calculate new height
  //for (int i=numNodes-1; i>=0; --i)
  float new_zeta;
  for (int i=0; i<numNodes; ++i)
  {

    // get the information about node relationships from the flow info object
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DR2;


    // check if this is a baselevel node
    if(node == receiver)
    {
      //cout << "This is a base level node. I don't update this node." << endl;
    }
    else if(elevation_map.find(node) != elevation_map.end())
    {
      //cout << "This is one of the fixed channel nodes. I don't update this node." << endl;
      FP_value = elevation_map[node];
      zeta[row][col]= FP_value;
    }
    else
    {
      // Get the data from the individual points
      K_value = K_values.get_data_element(row,col);
      U_value = U_values.get_data_element(row,col);

      // get the distance between nodes. Depends on flow direction
      switch (flow.retrieve_flow_length_code_of_node(node))
      {
        case 0:
          dx = -99;
        break;
        case 1:
          dx = DataResolution;
        break;
        case 2:
          dx = dx_root2;
        break;
        default:
          dx = -99;
        break;
      }


      receiver_elev = zeta[receiver_row][receiver_col];
      area_pow = pow(drainageArea,m_exp);
      parenth_term = U_value/(K_value*area_pow);
      new_zeta = receiver_elev+ dx*( pow(parenth_term,one_over_n));

      // now check to make sure the new zeta is not above the fixed channel


    }
    



  }
  //return LSDRasterModel(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This takes the model and calculates the steady state fluvial surface derived from
// a chi map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
float LSDRasterModel::fluvial_snap_to_steady_state_tune_K_for_relief(float U, float desired_relief)
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // We need to fill the raster so we don't get internally drained catchments
  float slope_for_fill = 0.0001;
  LSDRaster filled = temp.fill(slope_for_fill);
  LSDFlowInfo flow(boundary_conditions, filled);

  float K = get_K();
  float m_exp = get_m();
  float n_exp = get_n();
  float area_threshold = 0;
  float m_over_n = m_exp/n_exp;
  float A_0 = 1;
  float one_over_n = 1/n_exp;

  LSDRaster ChiValues = flow.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0,area_threshold);
  float thisChi;
  float MaxChi = 0;
  // now we calculate the elevations assuming that the elevation at chi = 0 is 0
  // This is based on equation 4 from mudd et al 2014 JGR-ES
  for (int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      thisChi = ChiValues.get_data_element(row,col);
      if (thisChi != NoDataValue)
      {
        if (thisChi > MaxChi)
        {
          MaxChi = thisChi;
        }
      }
    }
  }

  // calculate K (you need to do some algebra on equation 4 from mudd et al 2014 JGR-ES)
  K = U * pow((desired_relief/MaxChi),-n_exp);
  cout << "I am updating K to " << K << " to get your desired relief." << endl;


  // now recalculate the elevations
  // This is based on equation 4 from mudd et al 2014 JGR-ES
  for (int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      thisChi = ChiValues.get_data_element(row,col);
      if (thisChi == NoDataValue)
      {
        zeta[row][col] = NoDataValue;
      }
      else
      {
        if (n_exp == 1)
        {
          zeta[row][col] = (U/K)*thisChi;
        }
        else
        {
          zeta[row][col] = pow( (U/K),one_over_n)*thisChi;
        }

      }
    }
  }

  set_K(K);
  RasterData = zeta.copy();
  return(K);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This takes the model and calculates the K needed for a fixed relief
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
float LSDRasterModel::fluvial_calculate_K_for_steady_state_relief(float U, float desired_relief)
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // We need to fill the raster so we don't get internally drained catchments
  float slope_for_fill = 0.0001;
  LSDRaster filled = temp.fill(slope_for_fill);
  LSDFlowInfo flow(boundary_conditions, filled);

  float K = get_K();
  float m_exp = get_m();
  float n_exp = get_n();
  float area_threshold = 0;
  float m_over_n = m_exp/n_exp;
  float A_0 = 1;

  LSDRaster ChiValues = flow.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0,area_threshold);
  float thisChi;
  float MaxChi = 0;
  // now we calculate the elevations assuming that the elevation at chi = 0 is 0
  // This is based on equation 4 from mudd et al 2014 JGR-ES
  for (int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      thisChi = ChiValues.get_data_element(row,col);
      if (thisChi != NoDataValue)
      {
        if (thisChi > MaxChi)
        {
          MaxChi = thisChi;
        }
      }
    }
  }

  // calculate K (you need to do some algebra on equation 4 from mudd et al 2014 JGR-ES)
  K = U * pow((desired_relief/MaxChi),-n_exp);
  cout << "I am updating K to " << K << " to get your desired relief." << endl;

  return(K);
}





//============================================================================================================================
// WORKING HERE
// TRYING TO GET CRITICAL SLOPES WORKING
//============================================================================================================================


// Just a structure that define a node by ID and elevation
// Same principle that Martin's filling algorithm
struct MyNode
{
  float elevation;
  std::pair<int,int> ID;
};
bool operator>( const MyNode& lhs, const MyNode& rhs )
{
  return lhs.elevation > rhs.elevation;
}
bool operator<( const MyNode& lhs, const MyNode& rhs )
{
  return lhs.elevation < rhs.elevation;
}



LSDRaster LSDRasterModel::basic_valley_fill_critical_slope(float critical_slope, int contributing_pixel_threshold)
{
  
  Array2D<float> zeta=RasterData.copy();
  float sqrt2 = sqrt(2);

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);  

  // We need to fill the raster so we don't get internally drained catchments
  cout << "I'm getting the flow info" << endl;
  float slope_for_fill = 0.0001;
  LSDRaster filled = temp.fill(slope_for_fill);
  LSDFlowInfo flow(boundary_conditions, filled);
  int contributing_pixels;

  vector<int> S_vec = flow.get_SVector();

  // also I want to keep track on which nodes have been processed
  map<pair<int,int> ,bool> has_been_proc;
  //And a map that check whether a node is in the river
  map<pair<int,int> ,bool> is_river;
  // This tests the exsitence of the node in a stack
  map<pair<int,int> ,bool> is_in_DEM;
  // Let's initialise my output as well
  Array2D<float> output(NRows,NCols, NoDataValue);
  // Before processing my nodes, I need to set up my priority queue
  // And set 2 priority queues that will work together
  priority_queue< MyNode, vector<MyNode>, greater<MyNode> > PriorityQueue_1, PriorityQueue_2;
  // I'll need a switch to know which PQ is getting filled
  int working_PQ = 1;
  
  
  // let's go through the stack and get the river nodes
  // I'll fill my PQ1 with the river_nodes
  cout << "Filling the priority queues. The no data value is: " << NoDataValue <<endl;
  //cout << "NRows: " << NRows << " and NCols: " << NCols << endl;
  int this_row,this_col, current_node;
  for(size_t i =0; i<S_vec.size(); i++)
  {
    // Row and col of the next element from the stack, starting from the baselevel of the given stack
    current_node = S_vec[i];
    flow.retrieve_current_row_and_col(current_node,this_row,this_col);
    
    // Topo and drainage area stack
    float this_elev =temp.get_data_element(this_row,this_col);
    contributing_pixels = flow.retrieve_contributing_pixels_of_node(current_node);
    pair<int,int> this_node = {this_row, this_col};

    is_in_DEM[this_node] = true;
    
    // Initialising my node
    MyNode this_MyNode;
    this_MyNode.ID = this_node;
    this_MyNode.elevation = this_elev;

    // is it a river
    if(contributing_pixels>=contributing_pixel_threshold)
    {
      PriorityQueue_1.push(this_MyNode);
      is_river[this_node] = true;
      output[this_row][this_col] = this_elev;
      has_been_proc[this_node] = true;
    }
    else
    {
      is_river[this_node] = false;
      output[this_row][this_col] = NoDataValue;
      has_been_proc[this_node] = false;
    }
  }
  // std::cout << "is empty??"
  
  //cout << "Now it is time to create the slopes" << endl;
  // Alright, I have my list of node ordered by elevation (thanks to fastscape and stuff)
  while(PriorityQueue_1.empty() == false || PriorityQueue_2.empty() == false)
  {
    //cout << "Whohoo, starting with the priority queueueueueue" << endl;
    // I'll want my node
    pair<int,int> this_node;
    MyNode this_MyNode;
    if(working_PQ == 1)
    {
      this_MyNode = PriorityQueue_1.top();
      PriorityQueue_1.pop();
    }
    else
    {
      this_MyNode = PriorityQueue_2.top();
      PriorityQueue_2.pop();
    }
    this_node = this_MyNode.ID;
    float this_elevation = this_MyNode.elevation;
    int this_row = this_node.first;
    int this_col = this_node.second;
    // let me go through the neighbors
    //cout << "z: " << this_elevation << " r: " << this_row << " c: " << this_col << endl;
    vector<int> row_neighbors = {this_row-1,this_row,this_row+1};
    vector<int> col_neighbors = {this_col-1,this_col,this_col+1};
    for(int nR=0;nR<int(row_neighbors.size());nR++)
    {
      for(int nC=0;nC<int(col_neighbors.size());nC++)
      {
        // ignore that node
        if(nR ==1 && nC == 1)
        {
          continue;
        }
        // getting this neighbor
        int neighb_row = row_neighbors[nR];
        int neighb_col = col_neighbors[nC];



        // Make sure it is not at the edge of the DEM
        if (neighb_row >= 0 && neighb_row <= NRows - 1 && neighb_col >= 0 && neighb_col <= NCols - 1)
        {

    
          //cout << "nr: " << neighb_row << " nc: " << neighb_col << endl;
          pair<int,int> nenode = {neighb_row,neighb_col};

          if(not is_in_DEM[nenode])
          {
            //cout << "Hey this neighbour isn't in the DEM" << endl;
          }  
          else
          {
            if(not is_river[nenode] && not has_been_proc[nenode])
            {
              // If we got this far it means the node needs to be processed. 
              float dx;
              if((nR == 0 || nR == 2) && (nC == 0 || nC == 2))
              {
                dx = sqrt2*DataResolution;
              }
              else
              {
                dx = DataResolution;
              }
              // backcalculating the slope
              //cout << "I'm getting a slope at r: " << neighb_row << " c: " << neighb_col << " elev: " << this_elevation << endl; 
              float new_elev = critical_slope * dx + this_elevation;
              //cout << "1..";
              output[neighb_row][neighb_col] = new_elev;
              //cout << "2..";
              has_been_proc[nenode] = true;
              //cout << "done" <<endl;
              MyNode nextnode; 
              nextnode.ID = nenode;
              nextnode.elevation = new_elev;
              if(working_PQ == 1)
              {
                PriorityQueue_2.push(nextnode);
                // std::cout << "didd" << std::endl;
              }
              else
              {
                PriorityQueue_1.push(nextnode);
              }
            }       
          }
        }
      }
      if(working_PQ == 1 && PriorityQueue_1.empty() == true)
      {
        working_PQ = 2;
      }
      else if(working_PQ == 2 && PriorityQueue_2.empty() == true)
      {
        working_PQ = 1;
      }
    }
  }
  // for(size_t i=0)

  cout << "Right, all finished. Generating the output raster" << endl;
  LSDRaster output_raster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, output, GeoReferencingStrings);  
  return output_raster;
}


LSDRaster LSDRasterModel::basic_smooth(float central_pixel_weighting)
{

  // at the moment the boundary type can only be 0 and this is a periodic
  // boundary type at the E and W boundaries.

  Array2D<float> new_data(NRows,NCols,NoDataValue);
  float total_weighting;
  float total_sum;
  int rp1, rm1,cp1, cm1;
  int boundary_type = 0;    // This later allows the code to be flexible in terms of the boundary type. 
                            // Currently only periodic boundaries are on offer. 


  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      total_weighting = 0;
      total_sum = 0;

      rp1 = row+1;
      rm1 = row-1;
      cp1 = col+1;
      cm1 = col-1;

      // implement boundary conditions.
      if(boundary_type == 0)
      {
        if (rp1 == NRows)
        {
          rp1 = rm1;
        }
        if (rm1 == -1)
        {
          rm1 = rp1;
        }
        if (cp1 == NCols)
        {
          cp1 = 0;
        }
        if(cm1 == -1)
        {
          cm1 = NCols-1;
        }
      }
      else
      {
        if (rp1 == NRows)
        {
          rp1 = rm1;
        }
        if (rm1 == -1)
        {
          rm1 = rp1;
        }
        if (cp1 == NCols)
        {
          cp1 = 0;
        }
        if(cm1 == -1)
        {
          cm1 = NCols-1;
        }
      }

      if( RasterData[row][col] != NoDataValue)
      {
        total_weighting += central_pixel_weighting;
        total_sum += central_pixel_weighting*RasterData[row][col];

        // now go through all the other directions.
        if (RasterData[row][cp1] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[row][cp1];
        }
        if (RasterData[row][cm1] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[row][cm1];
        }
        if (RasterData[rp1][col] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[rp1][col];
        }
        if (RasterData[rm1][col] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[rm1][col];
        }
      }
      // Now update the array

      new_data[row][col] = total_sum/total_weighting;
    }
  }

  
  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, new_data, GeoReferencingStrings);  
  return temp;
}



LSDRaster LSDRasterModel::basic_valley_fill_critical_slope(LSDRaster& S_c_raster, int contributing_pixel_threshold)
{
  
  Array2D<float> zeta=RasterData.copy();
  float sqrt2 = sqrt(2);

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);  

  // We need to fill the raster so we don't get internally drained catchments
  cout << "I'm getting the flow info" << endl;
  float slope_for_fill = 0.0001;
  LSDRaster filled = temp.fill(slope_for_fill);
  LSDFlowInfo flow(boundary_conditions, filled);
  int contributing_pixels;

  vector<int> S_vec = flow.get_SVector();

  // also I want to keep track on which nodes have been processed
  map<pair<int,int> ,bool> has_been_proc;
  //And a map that check whether a node is in the river
  map<pair<int,int> ,bool> is_river;
  // This tests the exsitence of the node in a stack
  map<pair<int,int> ,bool> is_in_DEM;
  // Let's initialise my output as well
  Array2D<float> output(NRows,NCols, NoDataValue);
  // Before processing my nodes, I need to set up my priority queue
  // And set 2 priority queues that will work together
  priority_queue< MyNode, vector<MyNode>, greater<MyNode> > PriorityQueue_1, PriorityQueue_2;
  // I'll need a switch to know which PQ is getting filled
  int working_PQ = 1;
  
  
  // let's go through the stack and get the river nodes
  // I'll fill my PQ1 with the river_nodes
  cout << "Filling the priority queues. The no data value is: " << NoDataValue <<endl;
  //cout << "NRows: " << NRows << " and NCols: " << NCols << endl;
  int this_row,this_col, current_node;
  for(size_t i =0; i<S_vec.size(); i++)
  {
    // Row and col of the next element from the stack, starting from the baselevel of the given stack
    current_node = S_vec[i];
    flow.retrieve_current_row_and_col(current_node,this_row,this_col);
    
    // Topo and drainage area stack
    float this_elev =temp.get_data_element(this_row,this_col);
    contributing_pixels = flow.retrieve_contributing_pixels_of_node(current_node);
    pair<int,int> this_node = {this_row, this_col};

    is_in_DEM[this_node] = true;
    
    // Initialising my node
    MyNode this_MyNode;
    this_MyNode.ID = this_node;
    this_MyNode.elevation = this_elev;

    // is it a river
    if(contributing_pixels>=contributing_pixel_threshold)
    {
      PriorityQueue_1.push(this_MyNode);
      is_river[this_node] = true;
      output[this_row][this_col] = this_elev;
      has_been_proc[this_node] = true;
    }
    else
    {
      is_river[this_node] = false;
      output[this_row][this_col] = NoDataValue;
      has_been_proc[this_node] = false;
    }
  }
  // std::cout << "is empty??"
  
  //cout << "Now it is time to create the slopes" << endl;
  // Alright, I have my list of node ordered by elevation (thanks to fastscape and stuff)
  while(PriorityQueue_1.empty() == false || PriorityQueue_2.empty() == false)
  {
    //cout << "Whohoo, starting with the priority queueueueueue" << endl;
    // I'll want my node
    pair<int,int> this_node;
    MyNode this_MyNode;
    if(working_PQ == 1)
    {
      this_MyNode = PriorityQueue_1.top();
      PriorityQueue_1.pop();
    }
    else
    {
      this_MyNode = PriorityQueue_2.top();
      PriorityQueue_2.pop();
    }
    this_node = this_MyNode.ID;
    float this_elevation = this_MyNode.elevation;
    int this_row = this_node.first;
    int this_col = this_node.second;
    // let me go through the neighbors
    //cout << "z: " << this_elevation << " r: " << this_row << " c: " << this_col << endl;
    vector<int> row_neighbors = {this_row-1,this_row,this_row+1};
    vector<int> col_neighbors = {this_col-1,this_col,this_col+1};
    for(int nR=0;nR<int(row_neighbors.size());nR++)
    {
      for(int nC=0;nC<int(col_neighbors.size());nC++)
      {
        // ignore that node
        if(nR ==1 && nC == 1)
        {
          continue;
        }
        // getting this neighbor
        int neighb_row = row_neighbors[nR];
        int neighb_col = col_neighbors[nC];



        // Make sure it is not at the edge of the DEM
        if (neighb_row >= 0 && neighb_row <= NRows - 1 && neighb_col >= 0 && neighb_col <= NCols - 1)
        {

    
          //cout << "nr: " << neighb_row << " nc: " << neighb_col << endl;
          pair<int,int> nenode = {neighb_row,neighb_col};

          if(not is_in_DEM[nenode])
          {
            //cout << "Hey this neighbour isn't in the DEM" << endl;
          }  
          else
          {
            if(not is_river[nenode] && not has_been_proc[nenode])
            {
              // If we got this far it means the node needs to be processed. 
              float dx;
              if((nR == 0 || nR == 2) && (nC == 0 || nC == 2))
              {
                dx = sqrt2*DataResolution;
              }
              else
              {
                dx = DataResolution;
              }
              // backcalculating the slope
              //cout << "I'm getting a slope at r: " << neighb_row << " c: " << neighb_col << " elev: " << this_elevation << endl; 
              float new_elev = S_c_raster.get_data_element(neighb_row,neighb_col) * dx + this_elevation;
              //cout << "1..";
              output[neighb_row][neighb_col] = new_elev;
              //cout << "2..";
              has_been_proc[nenode] = true;
              //cout << "done" <<endl;
              MyNode nextnode; 
              nextnode.ID = nenode;
              nextnode.elevation = new_elev;
              if(working_PQ == 1)
              {
                PriorityQueue_2.push(nextnode);
                // std::cout << "didd" << std::endl;
              }
              else
              {
                PriorityQueue_1.push(nextnode);
              }
            }       
          }
        }
      }
      if(working_PQ == 1 && PriorityQueue_1.empty() == true)
      {
        working_PQ = 2;
      }
      else if(working_PQ == 2 && PriorityQueue_2.empty() == true)
      {
        working_PQ = 1;
      }
    }
  }
  // for(size_t i=0)

  cout << "Right, all finished. Generating the output raster" << endl;
  LSDRaster output_raster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, output, GeoReferencingStrings);  
  return output_raster;
}












//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates a new elevation raster based on instantaneous tilting of the model
// by a defined tilt angle.
// Need to specify: angle = tilt angle
//                  tilt_boundary = which boundary is the model tilted from. N, S, E, or W
// FJC 11/03/19
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::instantaneous_tilt(float angle, string tilt_boundary)
{
  Array2D<float> elev = RasterData.copy();
  Array2D<float> new_elev(NRows, NCols, NoDataValue);

  // loop through each possible tilt direction and get the new array of elevations
  // after tilting
  if (tilt_boundary == "N")  // north is tilt boundary - max elevation at the S
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // first row stays at same elevation
        if (i == 0) { new_elev[i][j] = elev[i][j]; }
        // other rows, calculate based on angle
        else
        {
          float this_elev = elev[i][j];
          float length = (i + 1) * DataResolution;
          new_elev[i][j] = (length * tan(angle)) + this_elev;
          cout << "old elev: " << this_elev << " new elev: " << (length * tan(angle)) + this_elev << endl;
        }
      }
    }
  }
  else if (tilt_boundary == "S")  // south is tilt boundary - max elevation at the N
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // last row stays at same elevation
        if (i == NRows-1) { new_elev[i][j] = elev[i][j]; }
        // other rows, calculate based on angle
        else
        {
          float this_elev = elev[i][j];
          float length = (NRows - i) * DataResolution;
          new_elev[i][j] = (length * tan(angle)) + this_elev;
          cout << "old elev: " << this_elev << " new elev: " << (length * tan(angle)) + this_elev << endl;
        }
      }
    }
  }
  else if (tilt_boundary == "E")  // east is tilt boundary - max elevation at the W
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // last col stays at same elevation
        if (j == NCols-1) { new_elev[i][j] = elev[i][j]; }
        // other cols, calculate based on angle
        else
        {
          float this_elev = elev[i][j];
          float length = (NCols - i) * DataResolution;
          new_elev[i][j] = (length * tan(angle)) + this_elev;
          cout << "old elev: " << this_elev << " new elev: " << (length * tan(angle)) + this_elev << endl;
        }
      }
    }
  }
  else if (tilt_boundary == "W")  // west is tilt boundary - max elevation at the E
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // first col stays at same elevation
        if (j == 0) { new_elev[i][j] = elev[i][j]; }
        // other cols, calculate based on angle
        else
        {
          float this_elev = elev[i][j];
          float length = (j + 1) * DataResolution;
          new_elev[i][j] = (length * tan(angle)) + this_elev;
          cout << "old elev: " << this_elev << " new elev: " << (length * tan(angle)) + this_elev << endl;
        }
      }
    }
  }
  else
  {
    cout << "Warning - you haven't set your boundary to N, W, E, or S. Returning the original raster" << endl;
    new_elev = elev;
  }

  // set the model to the array of new elevations
  RasterData = new_elev;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This is the component of the model that is solved using the
// FASTSCAPE algorithm of Willett and Braun (2013)
// Uses Newton's method to solve incision if the slope exponent != 1
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::fluvial_incision( void )
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  LSDFlowInfo flow(boundary_conditions, temp);
  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx, streamPowerFactor;
  float K = get_K();

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;


  // this is only for bug checking
  if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
  {
    cout << "Drainage area: " << endl;
    for (int i=0; i<NRows*NCols; ++i)
    {
      drainageArea = flow.retrieve_contributing_pixels_of_node(i) *  DR2;
      cout << drainageArea << " ";
      if (((i+1)%NCols) == 0)
      cout << endl;
    }
  }

  // Step two calculate new height
  //for (int i=numNodes-1; i>=0; --i)
  for (int i=0; i<numNodes; ++i)
  {

    // get the information about node relationships from the flow info object
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DR2;

    // some code for debugging
    if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
    {
      cout << row << ", " << col << ", " << receiver_row << ", " << receiver_col << endl;
      cout << flow.retrieve_flow_length_code_of_node(node) << endl;
      cout << drainageArea << endl;
    }

    // get the distance between nodes. Depends on flow direction
    switch (flow.retrieve_flow_length_code_of_node(node))
    {
      case 0:
  dx = -99;
  break;
      case 1:
        dx = DataResolution;
  break;
      case 2:
        dx = dx_root2;
  break;
      default:
  dx = -99;
  break;
    }

    // some logic if n is close to 1. Saves a bit of computational expense.
    if (abs(n - 1) < 0.0001)
    {
      if (dx == -99)
        continue;

      // compute new elevation if node is not a base level node
      if (node != receiver)
      {
        streamPowerFactor = K * pow(drainageArea, m) * (timeStep / dx);
        zeta[row][col] = (zeta[row][col] + zeta[receiver_row][receiver_col] * streamPowerFactor) /
                         (1 + streamPowerFactor);

        // check for overexcavation
        if(zeta[row][col] < zeta[receiver_row][receiver_col])
        {
          //cout << "Warning, overexcavation. Setting to minimum slope." << endl;
          zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
        }
      }
    }
    else    // this else loop is for when n is not close to one and you need an iterative solution
    {
      if (dx == -99)
        continue;
      float new_zeta = zeta[row][col];
      float old_zeta = zeta[row][col];

      float epsilon;     // in newton's method, z_n+1 = z_n - f(z_n)/f'(z_n)
                         // and here epsilon =   f(z_n)/f'(z_n)
                         // f(z_n) = -z_n + z_old - dt*K*A^m*( (z_n-z_r)/dx )^n
                         // We differentiate the above equation to get f'(z_n)
                         // the resulting equation f(z_n)/f'(z_n) is seen below
      float streamPowerFactor = K * pow(drainageArea, m) * timeStep;
      float slope;

      // iterate until you converge on a solution. Uses Newton's method.
      int iter_count = 0;
      do
      {
        slope = (new_zeta - zeta[receiver_row][receiver_col]) / dx;

        if(slope < 0)
        {
          epsilon = 0;
        }
        else
        {
          epsilon = (new_zeta - old_zeta + streamPowerFactor * pow(slope, n)) /
               (1 + streamPowerFactor * (n/dx) * pow(slope, n-1));
        }
        new_zeta -= epsilon;

        // This limits the number of iterations
        iter_count++;
        if(iter_count > 100)
        {
          //cout << "Too many iterations! epsilon is: " << abs(epsilon) << endl;
          epsilon = 0.5e-6;
        }
      } while (abs(epsilon) > 1e-6);
      zeta[row][col] = new_zeta;

      // check for overexcavation
      if(zeta[row][col] < zeta[receiver_row][receiver_col])
      {
        //cout << "Warning, overexcavation. Setting to minimum slope." << endl;
        zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
      }

    }
  }
  //return LSDRasterModel(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  this->RasterData = zeta.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This is the component of the model that is solved using the
// FASTSCAPE algorithm of Willett and Braun (2013)
// Uses Newton's method to solve incision if the slope exponent != 1
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::fluvial_incision_with_uplift( void )
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  //cout << "The nodatavalue is: " << NoDataValue << endl;
  LSDFlowInfo flow(boundary_conditions, temp);

  //for(int i = 0; i<4; i++)
  //{
  //  cout << "bc["<<i<<"]: " << boundary_conditions[i] << endl;
  //}

  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx, streamPowerFactor;
  float K = get_K();
  float U;

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;


  // this is only for bug checking
  if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
  {
    cout << "Drainage area: " << endl;
    for (int i=0; i<NRows*NCols; ++i)
    {
      drainageArea = flow.retrieve_contributing_pixels_of_node(i) *  DR2;
      cout << drainageArea << " ";
      if (((i+1)%NCols) == 0)
      cout << endl;
    }
  }

  // Step two calculate new height
  //for (int i=numNodes-1; i>=0; --i)
  for (int i=0; i<numNodes; ++i)
  {

    // get the information about node relashionships from the flow info object
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DR2;

    // some code for debugging
    if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
    {
      cout << row << ", " << col << ", " << receiver_row << ", " << receiver_col << endl;
      cout << flow.retrieve_flow_length_code_of_node(node) << endl;
      cout << drainageArea << endl;
    }

    // get the distance between nodes. Depends on flow direction
    switch (flow.retrieve_flow_length_code_of_node(node))
    {
      case 0:
        dx = -99;
        break;
      case 1:
        dx = DataResolution;
        break;
      case 2:
        dx = dx_root2;
        break;
      default:
        dx = -99;
        break;
    }

    // some logic if n is close to 1. Saves a bit of computational expense.
    if (abs(n - 1) < 0.0001)
    {
      if (dx == -99)
        continue;

      // compute new elevation if node is not a base level node
      if (node != receiver)
      {
        // get the uplift rate
        U = get_uplift_rate_at_cell(row,col);

        // get the stream power factor
        streamPowerFactor = K * pow(drainageArea, m) * (timeStep / dx);

        // calculate elevation
        zeta[row][col] = (zeta[row][col]
                          + zeta[receiver_row][receiver_col]*streamPowerFactor
                          + timeStep*U) /
                         (1 + streamPowerFactor);

        if(zeta[row][col] < zeta[receiver_row][receiver_col])
        {
          //cout << "Warning, overexcavation. Setting to minimum slope." << endl;
          zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
        }
      }
    }
    else    // this else loop is for when n is not close to one and you need an iterative solution
    {
      if (dx == -99)
      {
        //cout << "WTF, I am getting an invalid flow length code. LSDRastermodel 4523" << endl;
        continue;
      }
      float new_zeta = zeta[row][col];
      //float old_iter_zeta = zeta[row][col];
      float old_zeta = zeta[row][col];

      //cout << "computing for n != 1. " << endl;


      // get the uplift rate
      U = get_uplift_rate_at_cell(row,col);

      float epsilon;     // in newton's method, z_n+1 = z_n - f(z_n)/f'(z_n)
                         // and here epsilon =   f(z_n)/f'(z_n)
                         // f(z_n) = -z_n + z_old - dt*K*A^m*( (z_n-z_r)/dx )^n
                         // We differentiate the above equation to get f'(z_n)
                         // the resulting equation f(z_n)/f'(z_n) is seen below
      float streamPowerFactor = K * pow(drainageArea, m) * timeStep;
      float slope;

      // iterate until you converge on a solution. Uses Newton's method.
      int iter_count = 0;
      do
      {
        slope = (new_zeta - zeta[receiver_row][receiver_col]) / dx;

        if(slope < 0)
        {
          epsilon = 0;
        }
        else
        {
          // Get epsilon based on f(z_n)/f'(z_n)
          epsilon = (new_zeta - old_zeta
                     + streamPowerFactor * pow(slope, n) - timeStep*U) /
               (1 + streamPowerFactor * (n/dx) * pow(slope, n-1));
          //cout << "slope: " << slope << " epsilon: " << epsilon << endl;
        }

        new_zeta -= epsilon;

        iter_count++;
        if(iter_count > 100)
        {
          //cout << "Too many iterations! epsilon is: " << abs(epsilon) << endl;
          epsilon = 0.5e-6;
        }

      } while (abs(epsilon) > 1e-6);
      zeta[row][col] = new_zeta;

      // check for overexcavation
      if(zeta[row][col] < zeta[receiver_row][receiver_col])
      {
        //cout << "Warning, overexcavation. Setting to minimum slope." << endl;
        zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
      }
    }
  }

  //return LSDRasterModel(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  this->RasterData = zeta.copy();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This is the component of the model that is solved using the
// FASTSCAPE algorithm of Willett and Braun (2013)
// Uses Newton's method to solve incision if the slope exponent != 1
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterModel::fluvial_incision_with_uplift_and_variable_K( LSDRaster& K_raster )
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  //cout << "The nodatavalue is: " << NoDataValue << endl;
  LSDFlowInfo flow(boundary_conditions, temp);

  //for(int i = 0; i<4; i++)
  //{
  //  cout << "bc["<<i<<"]: " << boundary_conditions[i] << endl;
  //}

  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx, streamPowerFactor;
  float U;

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;


  // this is only for bug checking
  if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
  {
    cout << "Drainage area: " << endl;
    for (int i=0; i<NRows*NCols; ++i)
    {
      drainageArea = flow.retrieve_contributing_pixels_of_node(i) *  DR2;
      cout << drainageArea << " ";
      if (((i+1)%NCols) == 0)
      cout << endl;
    }
  }

  // Step two calculate new height
  //for (int i=numNodes-1; i>=0; --i)
  for (int i=0; i<numNodes; ++i)
  {

    // get the information about node relashionships from the flow info object
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DR2;

    // some code for debugging
    if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
    {
      cout << row << ", " << col << ", " << receiver_row << ", " << receiver_col << endl;
      cout << flow.retrieve_flow_length_code_of_node(node) << endl;
      cout << drainageArea << endl;
    }

    // get the distance between nodes. Depends on flow direction
    switch (flow.retrieve_flow_length_code_of_node(node))
    {
      case 0:
        dx = -99;
        break;
      case 1:
        dx = DataResolution;
        break;
      case 2:
        dx = dx_root2;
        break;
      default:
        dx = -99;
        break;
    }

    // some logic if n is close to 1. Saves a bit of computational expense.
    if (abs(n - 1) < 0.0001)
    {
      if (dx == -99)
        continue;

      // compute new elevation if node is not a base level node
      if (node != receiver)
      {
        // get the uplift rate
        U = get_uplift_rate_at_cell(row,col);

        // get the stream power factor
        streamPowerFactor = K_raster.get_data_element(row,col) * pow(drainageArea, m) * (timeStep / dx);

        // calculate elevation
        zeta[row][col] = (zeta[row][col]
                          + zeta[receiver_row][receiver_col]*streamPowerFactor
                          + timeStep*U) /
                         (1 + streamPowerFactor);

        if(zeta[row][col] < zeta[receiver_row][receiver_col])
        {
          zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
        }
      }
    }
    else    // this else loop is for when n is not close to one and you need an iterative solution
    {
      if (dx == -99)
      {
        continue;
      }
      float new_zeta = zeta[row][col];
      //float old_iter_zeta = zeta[row][col];
      float old_zeta = zeta[row][col];

      // get the uplift rate
      U = get_uplift_rate_at_cell(row,col);

      float epsilon;     // in newton's method, z_n+1 = z_n - f(z_n)/f'(z_n)
                         // and here epsilon =   f(z_n)/f'(z_n)
                         // f(z_n) = -z_n + z_old - dt*K*A^m*( (z_n-z_r)/dx )^n
                         // We differentiate the above equation to get f'(z_n)
                         // the resulting equation f(z_n)/f'(z_n) is seen below
      float streamPowerFactor = K_raster.get_data_element(row,col) * pow(drainageArea, m) * timeStep;
      float slope;

      // iterate until you converge on a solution. Uses Newton's method.
      int iter_count = 0;
      do
      {
        slope = (new_zeta - zeta[receiver_row][receiver_col]) / dx;

        if(slope < 0)
        {
          epsilon = 0;
        }
        else
        {
          // Get epsilon based on f(z_n)/f'(z_n)
          epsilon = (new_zeta - old_zeta
                     + streamPowerFactor * pow(slope, n) - timeStep*U) /
               (1 + streamPowerFactor * (n/dx) * pow(slope, n-1));
        }

        new_zeta -= epsilon;

        iter_count++;
        if(iter_count > 100)
        {
          epsilon = 0.5e-6;
        }

      } while (abs(epsilon) > 1e-6);
      zeta[row][col] = new_zeta;

      // check for overexcavation
      if(zeta[row][col] < zeta[receiver_row][receiver_col])
      {
        zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
      }
    }
  }
  this->RasterData = zeta.copy();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is the component of the model that is solved using the
// FASTSCAPE algorithm of Willett and Braun (2013)
// Uses Newton's method to solve incision if the slope exponent != 1
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::fluvial_incision_with_variable_uplift_and_variable_K( LSDRaster& Urate_raster, LSDRaster& K_raster )
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  //cout << "The nodatavalue is: " << NoDataValue << endl;
  LSDFlowInfo flow(boundary_conditions, temp);

  //for(int i = 0; i<4; i++)
  //{
  //  cout << "bc["<<i<<"]: " << boundary_conditions[i] << endl;
  //}

  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx, streamPowerFactor;
  float U;

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;


  // this is only for bug checking
  if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
  {
    cout << "Drainage area: " << endl;
    for (int i=0; i<NRows*NCols; ++i)
    {
      drainageArea = flow.retrieve_contributing_pixels_of_node(i) *  DR2;
      cout << drainageArea << " ";
      if (((i+1)%NCols) == 0)
      cout << endl;
    }
  }

  // Step two calculate new height
  //for (int i=numNodes-1; i>=0; --i)
  for (int i=0; i<numNodes; ++i)
  {

    // get the information about node relashionships from the flow info object
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DR2;

    // some code for debugging
    if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
    {
      cout << row << ", " << col << ", " << receiver_row << ", " << receiver_col << endl;
      cout << flow.retrieve_flow_length_code_of_node(node) << endl;
      cout << drainageArea << endl;
    }

    // get the distance between nodes. Depends on flow direction
    switch (flow.retrieve_flow_length_code_of_node(node))
    {
      case 0:
        dx = -99;
        break;
      case 1:
        dx = DataResolution;
        break;
      case 2:
        dx = dx_root2;
        break;
      default:
        dx = -99;
        break;
    }

    // some logic if n is close to 1. Saves a bit of computational expense.
    if (abs(n - 1) < 0.0001)
    {
      if (dx == -99)
        continue;

      // compute new elevation if node is not a base level node
      if (node != receiver)
      {
        // get the uplift rate
        U = Urate_raster.get_data_element(row,col);

        // get the stream power factor
        streamPowerFactor = K_raster.get_data_element(row,col) * pow(drainageArea, m) * (timeStep / dx);

        // calculate elevation
        zeta[row][col] = (zeta[row][col]
                          + zeta[receiver_row][receiver_col]*streamPowerFactor
                          + timeStep*U) /
                         (1 + streamPowerFactor);

        if(zeta[row][col] < zeta[receiver_row][receiver_col])
        {
          zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
        }
      }
    }
    else    // this else loop is for when n is not close to one and you need an iterative solution
    {
      if (dx == -99)
      {
        continue;
      }
      float new_zeta = zeta[row][col];
      //float old_iter_zeta = zeta[row][col];
      float old_zeta = zeta[row][col];

      // get the uplift rate
      U = Urate_raster.get_data_element(row,col);

      float epsilon;     // in newton's method, z_n+1 = z_n - f(z_n)/f'(z_n)
                         // and here epsilon =   f(z_n)/f'(z_n)
                         // f(z_n) = -z_n + z_old - dt*K*A^m*( (z_n-z_r)/dx )^n
                         // We differentiate the above equation to get f'(z_n)
                         // the resulting equation f(z_n)/f'(z_n) is seen below
      float streamPowerFactor = K_raster.get_data_element(row,col) * pow(drainageArea, m) * timeStep;
      float slope;

      // iterate until you converge on a solution. Uses Newton's method.
      int iter_count = 0;
      do
      {
        slope = (new_zeta - zeta[receiver_row][receiver_col]) / dx;

        if(slope < 0)
        {
          epsilon = 0;
        }
        else
        {
          // Get epsilon based on f(z_n)/f'(z_n)
          epsilon = (new_zeta - old_zeta
                     + streamPowerFactor * pow(slope, n) - timeStep*U) /
               (1 + streamPowerFactor * (n/dx) * pow(slope, n-1));
        }

        new_zeta -= epsilon;

        iter_count++;
        if(iter_count > 100)
        {
          epsilon = 0.5e-6;
        }

      } while (abs(epsilon) > 1e-6);
      zeta[row][col] = new_zeta;

      // check for overexcavation
      if(zeta[row][col] < zeta[receiver_row][receiver_col])
      {
        zeta[row][col] = zeta[receiver_row][receiver_col]+(0.00001)*dx;
      }
    }
  }
  this->RasterData = zeta.copy();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is the component of the model that is solved using the
// FASTSCAPE algorithm of Willett and Braun (2013)
// Uses Newton's method to solve incision if the slope exponent != 1
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::fluvial_incision_with_variable_uplift_and_variable_K_adaptive_timestep( LSDRaster& Urate_raster, LSDRaster& K_raster )
{
  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);
  //cout << "The nodatavalue is: " << NoDataValue << endl;
  LSDFlowInfo flow(boundary_conditions, temp);

  //for(int i = 0; i<4; i++)
  //{
  //  cout << "bc["<<i<<"]: " << boundary_conditions[i] << endl;
  //}

  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx, streamPowerFactor;
  float U;

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;

  // this is only for bug checking
  if (quiet == false && name == "debug" && NRows <= 10 && NCols <= 10)
  {
    cout << "Drainage area: " << endl;
    for (int i=0; i<NRows*NCols; ++i)
    {
      drainageArea = flow.retrieve_contributing_pixels_of_node(i) *  DR2;
      cout << drainageArea << " ";
      if (((i+1)%NCols) == 0)
      cout << endl;
    }
  }

  // We nest the main calculation routines in a logic statement that will
  // recalculate everything at a smaller timestep if the model overexcavates
  int timestep_iterator = 0;      // this checks how many times you have reduced the
                                  // timestep
  bool it_has_overexcavated;
  do
  {
    // reset the overexcavation switch
    it_has_overexcavated = false;

    // reset zeta to the old elevation. This wastes a bit of time but we need it
    // for the adaptive timestepping
    zeta=RasterData.copy();

    // Calculate new heights
    for (int i=0; i<numNodes; ++i)
    {
      // get the information about node relashionships from the flow info object
      node = nodeList[i];
      flow.retrieve_current_row_and_col(node, row, col);
      flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
      drainageArea = flow.retrieve_contributing_pixels_of_node(node)*DR2;

      // get the distance between nodes. Depends on flow direction
      switch (flow.retrieve_flow_length_code_of_node(node))
      {
        case 0:
          dx = -99;
          break;
        case 1:
          dx = DataResolution;
          break;
        case 2:
          dx = dx_root2;
          break;
        default:
          dx = -99;
          break;
      }

      // some logic if n is close to 1. Saves a bit of computational expense.
      if (abs(n - 1) < 0.0001)
      {
        if (dx == -99)
          continue;

        // compute new elevation if node is not a base level node
        if (node != receiver)
        {
          // get the uplift rate
          U = Urate_raster.get_data_element(row,col);

          // get the stream power factor
          streamPowerFactor = K_raster.get_data_element(row,col) * pow(drainageArea, m) * (timeStep / dx);

          // calculate elevation
          float zeta_old = zeta[row][col];
          zeta[row][col] = (zeta[row][col]
                          + zeta[receiver_row][receiver_col]*streamPowerFactor
                          + timeStep*U) /
                         (1 + streamPowerFactor);

          // check for overexcavation
          if(zeta[row][col] <= zeta[receiver_row][receiver_col])
          {
            //cout << "HEY HEY JABBA I found overexcavation!" << endl;
            it_has_overexcavated = true;

            if (timestep_iterator> 100)
            {
              cout << "There is an overexcavation that has not  been fixed by a very small timestep." << endl;
              cout << "I think there is a numerical instability and I am killing the computation." << endl;
              cout << "This does not mean that you are a bad person." << endl;
              cout << "zeta old is" << zeta_old << endl;
              cout << "zeta_ reciever is: " << zeta[receiver_row][receiver_col] << endl;
              cout << "Streampower factor: " << streamPowerFactor << endl;
              cout << "denominator is: " << (1 + streamPowerFactor) << endl;
              cout << "uplift term is: " << timeStep*U << endl;
              cout << "reciever term is: " << zeta[receiver_row][receiver_col]*streamPowerFactor << endl;
              exit(EXIT_FAILURE);
            }


          }
        }
      }
      else    // this else loop is for when n is not close to one and you need an iterative solution
      {
        if (dx == -99)
        {
          continue;
        }
        float new_zeta = zeta[row][col];
        //float old_iter_zeta = zeta[row][col];
        float old_zeta = zeta[row][col];

        // get the uplift rate
        U = Urate_raster.get_data_element(row,col);

        float epsilon;     // in newton's method, z_n+1 = z_n - f(z_n)/f'(z_n)
                          // and here epsilon =   f(z_n)/f'(z_n)
                            // f(z_n) = -z_n + z_old - dt*K*A^m*( (z_n-z_r)/dx )^n
                          // We differentiate the above equation to get f'(z_n)
                         // the resulting equation f(z_n)/f'(z_n) is seen below
        float streamPowerFactor = K_raster.get_data_element(row,col) * pow(drainageArea, m) * timeStep;
        float slope;

        // iterate until you converge on a solution. Uses Newton's method.
        int iter_count = 0;
        do
        {
          slope = (new_zeta - zeta[receiver_row][receiver_col]) / dx;

          if(slope < 0)
          {
            epsilon = 0;
          }
          else
          {
            // Get epsilon based on f(z_n)/f'(z_n)
            epsilon = (new_zeta - old_zeta
                     + streamPowerFactor * pow(slope, n) - timeStep*U) /
               (1 + streamPowerFactor * (n/dx) * pow(slope, n-1));
          }

          new_zeta -= epsilon;

          iter_count++;
          if(iter_count > 100)
          {
            epsilon = 0.5e-6;
          }

        } while (abs(epsilon) > 1e-6);
        zeta[row][col] = new_zeta;

        // check for overexcavation
        if(zeta[row][col] <= zeta[receiver_row][receiver_col])
        {
          it_has_overexcavated = true;

          // kill the program if the number of overexcavation steps get too small
          if (timestep_iterator> 100)
          {
            cout << "There is an overexcavation that has not  been fixed by a very small timestep." << endl;
            cout << "I think there is a numerical instability and I am killing the computation." << endl;
            cout << "This does not mean that you are a bad person." << endl;
            exit(EXIT_FAILURE);
          }
        }
      }        // end logic for n not equal to one

      if (it_has_overexcavated)
      {
        //cout << "Whoops I had an overexcavation! " << endl;
        //cout << " The number of times this has happend this timestep is: " << timestep_iterator << endl;
        timeStep = timeStep*0.25;
        timestep_iterator++;
        i = numNodes;
      }

    }         // end logic for node loop

  } while ( it_has_overexcavated);

  // if the model didn't overexcavate at all, then increase the timestep.
  //cout << "The timestep iterator is: " << timestep_iterator << endl;
  if(timestep_iterator ==0)
  {
    timeStep = 2*timeStep;
    if (timeStep > maxtimeStep)
    {
      timeStep = maxtimeStep;
    }
  }

  this->RasterData = zeta.copy();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is more or less identical to fluvial_incision above, but it
// Returns a raster with the erosion rate and takes arguments rather
// than reading from data members
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDRasterModel::fluvial_erosion_rate(float timestep, float K, float m, float n,
                                               vector <string> boundary)
{
  Array2D<float> erosionRate(NRows, NCols, NoDataValue);
  Array2D<float> zeta = RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDFlowInfo flow(boundary, *this);
  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx, streamPowerFactor;
  float root_2 = pow(2, 0.5);

  // Step two calculate new height
  for
   (int i=numNodes-1; i>=0; --i)
  {
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DataResolution * DataResolution;
    switch (flow.retrieve_flow_length_code_of_node(node))
    {
      case 0:
        dx = -99;
        break;
      case 1:
        dx = DataResolution;
        break;
      case 2:
        dx = DataResolution * root_2;
        break;
      default:
        dx = -99;
        break;
    }
    if (abs(n - 1) < 0.0001)
    {
      if (node == receiver)
      {
        erosionRate[row][col] = 0;
      }
      else
      {
      streamPowerFactor = K * pow(drainageArea, m) * (timestep / dx);
      erosionRate[row][col] = ((RasterData[row][col] + RasterData[receiver_row][receiver_col] *
              streamPowerFactor) /
             (1 + streamPowerFactor) - RasterData[row][col]) / timestep;
      }
    }
    else
    {
      float new_zeta = RasterData[row][col];
      float old_zeta = RasterData[row][col];

      float epsilon;
      float streamPowerFactor = K * pow(drainageArea, m) * timestep;
      float slope;
      do
      {
        slope = (new_zeta - zeta[receiver_row][receiver_col]) / dx;
        if (slope < 0)
        {
          epsilon = 0;
        }
        else
        {
          epsilon = (new_zeta - old_zeta + streamPowerFactor * pow(slope, n)) /
               (1 + streamPowerFactor * (n/dx) * pow(slope, n-1));
        }
        new_zeta -= epsilon;
      } while (abs(epsilon > 1e-6));
      erosionRate[row][col] = (new_zeta - old_zeta) / timestep;
    }
  }
  return LSDRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, erosionRate);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Returns the fluvial inicision rate, but uses data members
// also returns an array
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D<float> LSDRasterModel::fluvial_erosion_rate()
{
  Array2D<float> erosionRate(NRows, NCols, NoDataValue);
  Array2D<float> zeta = RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDFlowInfo flow(boundary_conditions, *this);
  vector <int> nodeList = flow.get_SVector();
  int numNodes = nodeList.size();
  int node, row, col, receiver, receiver_row, receiver_col;
  float drainageArea, dx, streamPowerFactor;
  float K = get_K();

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  float DR2 = DataResolution*DataResolution;

  // Step two calculate new height
  for (int i=numNodes-1; i>=0; --i)
  {
    node = nodeList[i];
    flow.retrieve_current_row_and_col(node, row, col);
    flow.retrieve_receiver_information(node, receiver, receiver_row, receiver_col);
    drainageArea = flow.retrieve_contributing_pixels_of_node(node) *  DR2;
    switch (flow.retrieve_flow_length_code_of_node(node))
    {
      case 0:
        dx = -99;
        break;
      case 1:
        dx = DataResolution;
        break;
      case 2:
        dx = dx_root2;
        break;
      default:
        dx = -99;
        break;
    }
    if (abs(n - 1) < 0.0001)
    {
      if (node == receiver)
      {
        erosionRate[row][col] = 0;
      }
      else
      {
      streamPowerFactor = K * pow(drainageArea, m) * (timeStep / dx);
      erosionRate[row][col] = ((RasterData[row][col] + RasterData[receiver_row][receiver_col] *
              streamPowerFactor) /
             (1 + streamPowerFactor) - RasterData[row][col]) / timeStep;
      }
    }
    else
    {
      float new_zeta = RasterData[row][col];
      float old_zeta = RasterData[row][col];

      float epsilon;
      float streamPowerFactor = K * pow(drainageArea, m) * timeStep;
      float slope;
      do
      {
        slope = (new_zeta - zeta[receiver_row][receiver_col]) / dx;
        epsilon = (new_zeta - old_zeta + streamPowerFactor * pow(slope, n)) /
             (1 + streamPowerFactor * (n/dx) * pow(slope, n-1));
        new_zeta -= epsilon;
      } while (abs(epsilon) > 1e-6);
      erosionRate[row][col] = (new_zeta - old_zeta) / timeStep;
    }
  }
  return erosionRate;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This little routine 'washes out' the channels: it assumes all
// sediment transported from hillslopes is transported  away in the channel.
// To determine the channel, it uses a threshold drainage area,
// which is one of the data members
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::wash_out( void )
{
  // don't do anything if there is no hillslope or fluvial erosion,
  // or if the threshold drainage is less than zero
  if (threshold_drainage < 0 || hillslope == false || fluvial == false)
    return;

  // get the old elevations
  Array2D<float> zeta=zeta_old.copy();
  int node;
  float DrainageArea;

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta);

  // SMM: does the flow info calculations every timestip: this might be streamlined
  // to do it every few timesteps
  LSDFlowInfo flow(boundary_conditions, temp);

  // loop through the nodes and wash out any node that is a channel
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      node = flow.retrieve_node_from_row_and_column(i, j);
      DrainageArea = flow.retrieve_contributing_pixels_of_node(node) *
                                               DataResolution * DataResolution;
      if (DrainageArea > threshold_drainage)
        RasterData[i][j] = zeta_old[i][j];
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




/*
LSDRasterModel LSDRasterModel::run_isostatic_correction( void )
{
  // Wrapper method for fortran program to calculate isostatic correction for
  // A given load on a plate
  // Original C program written by Jon Pelletier
  // Edited by James Jenkinson

  Array2D <float> output(NRows, NCols, NoDataValue);

  // Run Routine (Pelletier)
  flex2d(NRows, NCols, DataResolution, NoDataValue, RasterData, output);

  return LSDRasterModel(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, output);
}
*/





void LSDRasterModel::Airy_isostasy( void )
{
  // Density of materials (kg.m^{-3})
  float rho_c = 2650;
  float rho_m = 3300;

  float load;    // total load at each cell (elevation + root)
  float zeta_root;  // Height per depth of root
  zeta_root = (rho_m - rho_c) / rho_c;

  for (int i=0; i<NRows; ++i)
  {
    for (int j = 0; j<NCols; ++j)
    {
      load = RasterData[i][j] + root_depth[i][j];
      root_depth[i][j] = load / (1 + zeta_root);
      RasterData[i][j] = load - root_depth[i][j];
    }
  }
}

void LSDRasterModel::flexural_isostasy( float alpha )
{
  int iter=0, max_iter = 200;
  Array2D<float> old_root;
  Array2D<float> difference;
  float max_error;
  float epsilon=0.0001;
  stringstream ss;

  do {
    ++iter;
    max_error = 0;
    old_root = root_depth.copy();
    root_depth = calculate_root();
    difference = root_depth - old_root;

    if ( quiet == false && name == "debug" && NRows <= 10 && NCols<= 10)
    {
      cout << "Topography: " << endl;
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        cout << RasterData[i][j] << " ";
      }
      cout << endl;
    }
      cout << "Root: " << endl;
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        cout << root_depth[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
    cout << "Difference: " << endl;
    for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        cout << difference[i][j] << " ";
      }
      cout << endl;
    }
    }

    //RasterData += difference;
    for (int i = 0; i < NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
      {
        RasterData[i][j] -= (difference[i][j] * alpha);
        root_depth[i][j] = old_root[i][j] + (difference[i][j] * alpha);
        if (abs(difference[i][j]) > max_error)
          max_error = abs(difference[i][j]);
      }
    }
    //cout << max_error << endl;

    //if (name == "debug")
    //{
      ss.str("");
      ss << "step" << iter;
      write_root(ss.str(), "asc");
      ss.str("");
      ss << "step_raster" << iter;
      write_raster(ss.str(), "asc");
    //}

  } while (max_error > epsilon && iter < max_iter);
}

void LSDRasterModel::flexural_isostasy_alt( void )
{
  Array2D<float> old_root;
  Array2D<float> difference;

  old_root = root_depth.copy();
  root_depth = calculate_root();
  difference = root_depth - old_root;

  if (quiet == false && name == "debug" && NRows <= 10 && NCols<= 10)
  {
    cout << "Topography: " << endl;
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      cout << RasterData[i][j] << " ";
    }
    cout << endl;
  }
    cout << "Root: " << endl;
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      cout << root_depth[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
  cout << "Difference: " << endl;
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      cout << difference[i][j] << " ";
    }
    cout << endl;
  }
  }

  //RasterData += difference;
  for (int i = 0; i < NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      RasterData[i][j] -= (difference[i][j]);
      root_depth[i][j] = old_root[i][j] + (difference[i][j]);
    }
  }
}

Array2D <float> LSDRasterModel::calculate_root( void )
{
  // Calculate padding size
  int Ly = pow(2, ceil(log(NRows)/log(2)));
  int Lx = pow(2, ceil(log(NCols)/log(2)));

  // Array structures for fourier coefficients
  Array2D <float> real_coeffs(Ly, Lx);    // Real part of fourier coefficients
  Array2D <float> imag_coeffs(Ly, Lx);    // Imaginary part "                "
  Array2D <float> detrend(Ly, Lx, 0.0);  // Detrended raster data
  Array2D <float> trend(NRows, NCols);    // Trend plane
  Array2D <float> output(NRows, NCols, 0.0);  // Output Array

  // Detrend the data
  detrend2D(RasterData, detrend, trend);

  // Set up parameters
  /*
  float E = 10E9;    // young's modulus (kg.m^{-1}.s^{-2} || Pa
  float Te = 10E3;    // Elastic thickness (m)
  float v = 0.20;    // poisson ratio
  */
  //float D = E * pow(Te, 3) / (12 * (1 - v*v));
  float D = rigidity;

  float rho_c = 2650;    // density of crust (kg.m^3)
  float rho_m = 3300;    // density of mantle (kg.m^3)
  float pi = 3.14159265359;  // pi
  float g = 9.81;    // gravitational acceleration (m.s^{-2})

  //cout << D << endl;

  // Calculate fourier transfor coefficients of load
  dfftw2D_fwd(detrend, real_coeffs, imag_coeffs, -1);

  // Shift the dataset
  Array2D <float> real_shift(Ly, Lx);
  Array2D <float> imag_shift(Ly, Lx);
  shift_spectrum(real_coeffs, imag_coeffs, real_shift, imag_shift);

  // Multiply coefficients by function (see Peletier, 2008)
  float coeff;

  //cout << "  Solving isostasy with Vening-Meinesz method" << endl;
  for (int i=0; i<Ly; ++i)
  {
    for (int j=0; j<Lx; ++j)
    {
      //coeff = rho_m/rho_c - 1 + (D / (rho_c*g)) * pow(pi*(pow(((float)i/Ly)*((float)i/Ly)
      //     + ((float)j/Lx)*((float)j/Lx), 0.5)), 4);
      //coeff = (rho_c / (rho_m - rho_c)) / ((D * pow(pi*(pow(((float)i/Ly)*((float)i/Ly)
      //     + ((float)j/Lx)*((float)j/Lx), 0.5)), 4))/((rho_m-rho_c)*g) + 1);
      coeff = (rho_c/ (rho_m-rho_c)) / (1 + 4*(4*D/(pow((rho_m-rho_c)*g,0.5) * pow(pi*(pow(((float)i/Ly)*((float)i/Ly)
           + ((float)j/Lx)*((float)j/Lx), 0.5)), 4))));
      real_shift[i][j] *= coeff;
      imag_shift[i][j] *= coeff;
    }
  }

  // De-shift spectrum
  shift_spectrum_inv(real_shift, imag_shift, real_coeffs, imag_coeffs);

  // Reconstruct array from new fourier coefficients
  dfftw2D_inv(real_coeffs, imag_coeffs, detrend, 1);

  // Reapply trend plane
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      if (i==0 && boundary_conditions[0][0] == 'b')
        output[i][j] = 0;
      else if (j==0 && boundary_conditions[3][0] == 'b')
        output[i][j] = 0;
      else if (i==NRows-1 && boundary_conditions[2][0] == 'b')
        output[i][j] = 0;
      else if (j==NCols-1 && boundary_conditions[1][0] == 'b')
        output[i][j] = 0;
      else
        output[i][j] = detrend[i][j]/(Lx*Ly) + trend[i][j];
    }
  }

  // DEBUG: print out isostasy map
  //LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, output);
  //temp.write_raster("output", "asc");
  //root_depth = output;
  return output;
}

Array2D<float> LSDRasterModel::calculate_airy( void )
{
  float rho_c = 2650;    // density of crust (kg.m^3)
  float rho_m = 3300;    // density of mantle (kg.m^3)

  Array2D <float> airy(NRows, NCols, 0.0);

  for (int i = 0; i<NRows; ++i)
  {
    for (int j = 0; j <NCols; ++j)
    {
      airy[i][j] = RasterData[i][j] * rho_c/(rho_m-rho_c);
    }
  }
  return airy;
}

void LSDRasterModel::write_root( string name, string ext )
{
  LSDRaster root(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, root_depth);
  root.write_raster(name, ext);
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints parameters to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::print_parameters( void )
{
  if (quiet)
    return;
  cout << "\n========================================================" << endl;
  cout << "\nModel run: " << name << endl;
  cout << "\nFrom 0 to " << endTime << " years, in increments of " << timeStep << endl;
  cout << NRows << " by " << NCols << endl;
  cout << "Cells " << DataResolution << " metres wide." << endl;

  cout << "\n---------------------------------" << endl;
  cout << "Boundary conditions: " << endl;

  for (int i=0; i<4; ++i)
  {
    switch (i) {
      case 0:
        cout << "North:\t"; break;
      case 1:
        cout << "East:\t"; break;
      case 2:
        cout << "South:\t"; break;
      case 3:
        cout << "West:\t"; break;
    }
    if (boundary_conditions[i][0] == 'b')     cout << "Base level" << endl;
    else if (boundary_conditions[i][0] == 'p')   cout << "Periodic" << endl;
    else            cout << "No flow" << endl;
  }

  cout << "\n---------------------------------" << endl;
  if (fluvial)
  {
    cout << "Fluvial:\tOn" << endl;
    cout << "\nFLUVIAL PARAMETERS:" << endl;
    cout << "\tK:\t\t" << K_fluv << endl;
    cout << "\tm:\t\t" << m << endl;
    cout << "\tn:\t\t" << n << endl;
  }
  else
    cout << "\nFluvial:\tOff" << endl;

  cout << "\n---------------------------------" << endl;
  if (hillslope)
  {
    cout << "Hillslope:\tOn\t" << ((nonlinear) ? "Non-linear" : "Linear") << endl;
    cout << "\nSOIL PARAMTERS:" << endl;
    cout << "\tD:\t\t" << K_soil << endl;
    if (nonlinear)
      cout << "\tCritical slope:\t" << S_c << endl;
  }
  else
    cout << "Hillslope:\tOff" << endl;

  cout << "\n---------------------------------" << endl;
  cout << "\nIsostasy:\t" << ((isostasy) ? "On" : "Off") << endl;
  if (isostasy)
    cout << "\tModel:\t\t" << ((flexure) ? "Flexural" : "Airy") << endl;

  cout << "\n========================================================" << endl;
  cout << "\n" << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This method calculates some features of the landscape at set times. The
// frequency of the reports are set by the data member report_delay
// One of the things it does is calculates erosion rates, and stores this
// as a data member
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::write_report( void )
{
  static ofstream outfile;

  // check to see if enough time has elapsed to write the report
  if (reporting && current_time > report_delay)
  {
    if ( outfile.is_open() == false)
    {
      // Headers
      outfile.open((report_name + "_report").c_str());
      outfile << name << endl;
      outfile << "Time\t";
      outfile << "Periodicity\t";
      outfile << ((fluvial) ? "K\t" : "");
      outfile << ((hillslope) ? "D\t" : "");
      outfile << "Erosion\t";
      outfile << "Total erosion\t";
      outfile << "Steady\t";
      outfile << "Max_height\t";
      outfile << "Mean_height\t";
      outfile << "Relief-3px\t";
      outfile << "Relief-10m\t";
      //outfile << "Relief-30m\t";
      outfile << "Drainage-20m2\t";
      outfile << "Drainage-200m2\t";
      //outfile << "Drainage-500m2\t";
      outfile << endl;
    }
    if ( recording == false)
      check_recording();
    if (print_erosion_cycle)
      erosion_cycle_field = Array2D<float>(NRows, NCols, 0.0);

    outfile << current_time << "\t";
    outfile << periodicity << "\t";
    if (fluvial)   outfile << get_K() << "\t";
    if (hillslope)   outfile << get_D() << "\t";
  }

  // Calculate erosion across landscape
  erosion_last_step = erosion;     // reset the erosion rates
  erosion = 0.0;
  float e;
  int n = 0;
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      e = get_erosion_at_cell(i,j);
      if (print_erosion_cycle && ((initial_steady_state || cycle_steady_check) && (K_mode !=0 || D_mode != 0)))
      {
        erosion_cycle_field[i][j] += e;
      }
      if ( is_base_level(i,j) == false)
      {
        erosion += e;
        ++n;
      }
    }
  }
  //cout << "\nn: " << n << endl;
  erosion /= n;
  // Add to total erosion count
  if (recording)
    total_erosion += erosion;
  // Calculate local response
  if (erosion > erosion_last_step)
    max_erosion = erosion;
  else if (erosion < erosion_last_step)
    min_erosion = erosion;
  if (min_erosion != -99 && max_erosion - min_erosion > response)
    response = max_erosion - min_erosion;
  if (recording)
  {
    if (erosion > max_erosion)
      max_erosion = erosion;
    if (min_erosion == -99 || erosion < min_erosion)
      min_erosion = erosion;
  }
  float mean_elev, max_elev, relief0, relief10;
  max_elev = max_elevation();
  mean_elev = mean_elevation();
  relief0 = mean_relief(0);
  relief10 = mean_relief(10);
  if (reporting && current_time > report_delay)
  {
    outfile << erosion << "\t";
    outfile << total_erosion << "\t";
    outfile << steady_state << "\t";
    outfile << max_elev << "\t";
    outfile << mean_elev << "\t";
    outfile << relief0 << "\t";
    outfile << relief10 << "\t";
    //outfile << mean_relief(30) << "\t";
    //outfile << mean_drainageDensity(500);
    outfile << endl;
  }
  if ((initial_steady_state || cycle_steady_check) && (K_mode !=0 || D_mode != 0))
    cycle_report(mean_elev, relief0, relief10);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// SMM: This seems to be a method to check on what has happened over a complete
// cycle.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::cycle_report( float elev, float relief0, float relief10)
{
  // There's got to be a better way to design this method, I don't like it
  static ofstream outfile;
  static int phase_pos = 1;
  if ( outfile.is_open() == false && reporting && current_time > report_delay)
  {
    outfile.open((report_name + "_cycle_report").c_str());
    outfile << name << endl;
    outfile << "Cycle\t";
    outfile << "Start_time\t";
    outfile << "End_time\t";
    outfile << "Periodicity\t";
    outfile << "Erosion\t";
    outfile << "Erosion_response\t";
    outfile << "Elevation\t";
    outfile << "Elevation_response\t";
    outfile << "Relief-3px\t";
    outfile << "Relief-3px_response\t";
    outfile << "Relief-10m\t";
    outfile << "Relief-10m_response\t";
    outfile << "Drainage-20m2\t";
    outfile << "Drainage-20m2_response\t";
    outfile << "Drainage-200m2\t";
    outfile << "Drainage-200m2_response\t";
    outfile << endl;
  }
  static float mean_eros=0,  mean_elev=0,  mean_relief0=0,  mean_relief10=0;
  static float max_eros=0,   max_elev=0,   max_relief0=0,   max_relief10=0;
  static float min_eros=-99, min_elev=-99, min_relief0=-99, min_relief10=-99;
  static int n = 0;
  static float start_time = current_time;

  if (current_time == 0)
  {
      mean_eros=0,  mean_elev=0,  mean_relief0=0,  mean_relief10=0;
      max_eros=0,   max_elev=0,   max_relief0=0,   max_relief10=0;
      min_eros=-99, min_elev=-99, min_relief0=-99, min_relief10=-99;
      n= 0;
  }
  // Check next cycle number
  current_time += timeStep;
  float p = periodicity;

  // if we are checking if the method is checking for steady state over a cycle,
  // initiate the erosion_cycle_record data member
  if (cycle_steady_check)
  {
    if (erosion_cycle_record.size() == 0)
      erosion_cycle_record = vector<float> (5, -99);
  else
    erosion_cycle_record.empty();
  }

  // SMM this seems to use a number of statistics derived from the underlying
  // raster, but they are not calculated here. Will need to find where
  // they are calculated.
  if (periodic_parameter(1, 1) > 1)
  {
    if (phase_pos == 0)
    {
      ++cycle_number;
      if(reporting && current_time > report_delay)
      {
        outfile << cycle_number-1 << "\t";
        outfile << start_time << "\t";
        outfile << current_time-timeStep << "\t";
        outfile << p << "\t";
        outfile << mean_eros/n << "\t";
        outfile << max_eros-min_eros << "\t";
        outfile << mean_elev/n << "\t";
        outfile << max_elev - min_elev << "\t";
        outfile << mean_relief0/n << "\t";
        outfile << max_relief0-min_relief0 << "\t";
        outfile << mean_relief10/n << "\t";
        outfile << max_relief10-min_relief10 << "\t";
        outfile << endl;
        start_time = current_time-timeStep;
      }
      if (print_erosion_cycle)
      {
        for (int i=0; i<NRows; ++i)
        {
          for (int j=0; j<NCols; ++j)
          {
            erosion_cycle_field[i][j] /= n;
          }
        }
        stringstream ss;
        ss << name << cycle_number -1 << "_cycle_erosion";
        LSDRaster * e_cycle;
        e_cycle = new LSDRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, erosion_cycle_field);
        e_cycle->write_raster(ss.str(), "asc");
        delete e_cycle;
        erosion_cycle_field = Array2D<float>(NRows, NCols, 0.0);
      }

      if (cycle_steady_check)
      {
        for (int i=0; i<4; ++i)
          erosion_cycle_record[i] = erosion_cycle_record[i+1];
        erosion_cycle_record[4] = mean_eros/n;
      }

      mean_eros=0,  mean_elev=0,  mean_relief0=0,  mean_relief10=0;
      max_eros=0,   max_elev=0,   max_relief0=0,   max_relief10=0;
      min_eros=-99, min_elev=-99, min_relief0=-99, min_relief10=-99;
      n = 0;
    }
    phase_pos = 1;
  }
  else
  {
    phase_pos = 0;
  }
  current_time -= timeStep;
  mean_elev += elev;
  mean_eros += erosion;
  mean_relief0 += relief0;
  mean_relief10 += relief10;

  if (elev > max_elev) max_elev = elev;
  if (erosion > max_eros) max_eros = erosion;
  if (relief0 > max_relief0) max_relief0 = relief0;
  if (relief10 > max_relief10) max_relief10 = relief10;

  if (min_elev==-99 || elev<min_elev) min_elev=elev;
  if (min_eros==-99 || erosion<min_eros) min_eros=erosion;
  if (min_relief0==-99 || relief0<min_relief0) min_relief0=relief0;
  if (min_relief10==-99 || relief10<min_relief10) min_relief10=relief10;

  ++n;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//   Some overall details about runs through periods, used by JAJ's periodicity
//  routines
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::final_report( void )
{
  ofstream final_report;
  final_report.open((report_name + "_final").c_str());
  final_report << name << endl;

  float run_time;
  run_time = (K_mode != 0 || D_mode != 0) ? current_time - time_delay - periodicity : current_time - time_delay;
  //cout << "Run time: " << run_time << endl;


  final_report << "Erosion\tAveraged\tResponse\tK amp\tD amp\tPeriodicity\tOvershoot" << endl;
  final_report << total_erosion << "\t" << (total_erosion / ( run_time * num_runs)) << "\t"
    << ((initial_steady_state) ? response/num_runs : -99) << "\t" << K_amplitude << "\t"
    << D_amplitude << "\t" << periodicity << "\t" << current_time - endTime << endl;
  //final_report << "Erosion\tResponse\tK amp\tD amp\tPeriodicity" << endl;
  //final_report << total_erosion << "\t" << ((initial_steady_state) ? max_erosion - min_erosion : -99) << "\t" << K_amplitude << "\t"
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function writes some averaged information about the erosion and
// apparent erosion rates
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::print_average_erosion_and_apparent_erosion( int frame,
                                 vector<LSDParticleColumn>& CRNColumns,
                                 LSDCRNParameters& CRNParams)
{

  float avg_uplift = get_average_upflit_rate_last_timestep();
  float avg_erosion = get_total_erosion_rate_over_timestep();

  static ofstream er_outfile;
  // Print the cosmo metadata
  if (er_outfile.is_open() == false)
  {
    string metadata_fname =  name+".er_frame_metadata";
    cout << "Name of raster metadata file is: " <<  metadata_fname << endl;
          er_outfile.open(metadata_fname.c_str());
           er_outfile << name << endl;
          er_outfile << "Frame_num\t";
    er_outfile << "Time\t";
    er_outfile << "K\t";
    er_outfile << "D\t";
    er_outfile << "Erosion\t";
    er_outfile << "Avg_uplift\t";
    er_outfile << "10Be_apparent_erosion\t";
    er_outfile << "14C_apparent_erosion\t";
    er_outfile << "21Ne_apparent_erosion\t";
    er_outfile << endl;
  }

  // also print the cosmo properties
  string frame_name = itoa(frame);
  string uscore = "_";
  string cosmo_fname = "_cosmo.cdata";
  cosmo_fname = name+uscore+frame_name+cosmo_fname;
  ofstream cosmo_out;
  cosmo_out.open(cosmo_fname.c_str());

  float tot_weighted_erate_10Be = 0;
  float tot_weighted_erate_14C = 0;
  float tot_weighted_erate_21Ne = 0;
  float tot_erate = 0;

  int N_CRNcols = int(CRNColumns.size());
  for (int cc = 0; cc<N_CRNcols; cc++)
  {

    // get the actual erosion at the cell
    double this_erosion = get_erosion_at_cell(CRNColumns[cc].getRow(),
                                              CRNColumns[cc].getCol());

    // get the uplift rate
    double this_U = get_uplift_rate_at_cell(  CRNColumns[cc].getRow(),
                                              CRNColumns[cc].getCol());

    // get the diffusivity
    double this_D = get_D();

    // get the K value
    double this_K = get_K();

    // get the apparent erosion at the cell
    vector<double> this_app_erosion =
       CRNColumns[cc].calculate_app_erosion_3CRN_neutron_rock_only(CRNParams);

    // get the totals for weighted average
    tot_erate += this_erosion;
    tot_weighted_erate_10Be += this_erosion*this_app_erosion[0];
    tot_weighted_erate_14C += this_erosion*this_app_erosion[1];
    tot_weighted_erate_21Ne += this_erosion*this_app_erosion[2];


    // print these things
    cosmo_out << current_time << " "<< CRNColumns[cc].getRow()
              << " " << CRNColumns[cc].getCol()
              << " " << this_D << " " << this_K << " " << this_U
              << " " << this_erosion << " " << this_app_erosion[0] << " "
              << " " << this_app_erosion[1] << " " << this_app_erosion[2]
              << endl;



  }
  cosmo_out.close();

  float weight_10Be =  tot_weighted_erate_10Be/tot_erate;
  float weight_14C =  tot_weighted_erate_14C/tot_erate;
  float weight_21Ne =  tot_weighted_erate_21Ne/tot_erate;

  er_outfile << frame << "\t";
  er_outfile << current_time << "\t";
  er_outfile << get_K() << "\t";
  er_outfile << get_D() << "\t";
  er_outfile << avg_erosion << "\t";
  er_outfile << avg_uplift << "\t";
  er_outfile << weight_10Be << "\t";
  er_outfile << weight_14C << "\t";
  er_outfile << weight_21Ne << "\t";
  er_outfile <<  endl;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function writes all of the erosion data from individual columns
// Warning: Writes very big files!!
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::print_column_erosion_and_apparent_erosion( int frame,
                                 vector<LSDParticleColumn>& CRNColumns,
                                 LSDCRNParameters& CRNParams)
{

  float avg_uplift = get_average_upflit_rate_last_timestep();
  float avg_erosion = get_total_erosion_rate_over_timestep();

  // get the number of columns
  int N_CRNcols = int(CRNColumns.size());

  // This if statement opens the file if it doesn't exist
  static ofstream er_outfile;
  // Print the cosmo metadata
  if (er_outfile.is_open() == false)
  {
    string CRNdata_fname =  name+".CRN_frame_metadata";
    cout << "Name of CRN file is: " <<  CRNdata_fname << endl;
          er_outfile.open(CRNdata_fname.c_str());
           er_outfile << name << endl;
          er_outfile << "Frame_num\t";
    er_outfile << "Time\t";
    er_outfile << "K\t";
    er_outfile << "D\t";
    er_outfile << "Avg_uplift\t";
    er_outfile << "Avg_erosion\t";

    // Now you need to loop through every column collecting data
    for (int cc = 0; cc<N_CRNcols; cc++)
    {
      er_outfile << "Row\tCol\t\Erosion\t10BeErosion";
    }

    er_outfile << endl;
  }


  // if the file is aloready open, print the data to file
  er_outfile << frame << "\t";
  er_outfile << current_time << "\t";
  er_outfile << get_K() << "\t";
  er_outfile << get_D() << "\t";
  er_outfile << avg_uplift << "\t";
  er_outfile << avg_erosion << "\t";

  // loop through the columns, printing the local data
  for (int cc = 0; cc<N_CRNcols; cc++)
  {
    er_outfile << CRNColumns[cc].getRow() << "\t";
    er_outfile << CRNColumns[cc].getCol() << "\t";


    // get the actual erosion at the cell
    double this_erosion = get_erosion_at_cell(CRNColumns[cc].getRow(),
                                              CRNColumns[cc].getCol());

    // get the uplift rate
    //double this_U = get_uplift_rate_at_cell(  CRNColumns[cc].getRow(),
    //                                          CRNColumns[cc].getCol());

    // get the apparent erosion at the cell
    vector<double> this_app_erosion =
       CRNColumns[cc].calculate_app_erosion_3CRN_neutron_rock_only(CRNParams);

    er_outfile << this_erosion << "\t";
    er_outfile << this_app_erosion[0] << "\t";

  }
  er_outfile <<  endl;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This routine prints rasters according to some internal model switched
// Used by JAJ's code so DO NOT MODIFY
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::print_rasters( int frame )
{
  string outfile_format = "bil";

  cout << endl;
  //cout << "Printing raster, initial steady state: " << initial_steady_state
  //     << " and D mode is : " << D_mode << endl;
  //cout << "current time" << current_time << " td: " << time_delay << " sd: " << switch_delay << endl;
  //cout << "Time in wave: " << (current_time - time_delay - switch_delay) << endl;
  //cout << "Periodicity: " << periodicity << " and pi: " << PI << endl;
  //cout << "input to sine wave: " << (current_time - time_delay - switch_delay) * 2 * PI / periodicity << endl;
  //cout << "Sin wave:" << sin( (current_time - time_delay - switch_delay) * 2 * PI / periodicity )  << endl;
  static ofstream outfile;
  if ( outfile.is_open() == false)
  {
    string metadata_fname =  name+"._frame_metadata";
    cout << "Name of raster metadata file is: " <<  metadata_fname << endl;
          outfile.open(metadata_fname.c_str());
           outfile << name << endl;
          outfile << "Frame_num\t";
    outfile << "Time\t";
    outfile << "K\t";
    outfile << "D\t";
    outfile << "Erosion\t";
    outfile << "Max_uplift\t";
    outfile << endl;
  }
  outfile << frame << "\t";
  outfile << current_time << "\t";
  outfile << get_K() << "\t";
  outfile << get_D() << "\t";
  outfile << erosion << "\t";
  outfile << get_max_uplift(); // << "\t"; i think this make the output buggy, at least for me it sxrew up the csv file! - BG
  outfile << endl;

  map<string,string> GRS = get_GeoReferencingStrings();

  //cout << "Printing, print elevation is " << print_elevation
  //     << " and erosion is " << print_erosion << endl;

  stringstream ss;
  if (print_elevation)
  {
    ss << name << frame;
    this->write_raster(ss.str(), outfile_format);
  }
  if (print_hillshade)
  {
    cout << "Printing the hillshade" << endl;
    ss.str("");
    ss << name << frame << "_hs";
    LSDRaster * hillshade;
    hillshade = new LSDRaster(*this);
    *hillshade = this->hillshade(45, 315, 1);
    hillshade->write_raster(ss.str(), outfile_format);
    delete hillshade;
  }
  if (print_erosion)
  {
    ss.str("");
    ss << name << frame << "_erosion";
    Array2D <float> erosion_field = calculate_erosion_rates( );
    LSDRaster * erosion;
    erosion = new LSDRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, erosion_field,GRS);
    erosion->write_raster(ss.str(), outfile_format);
    delete erosion;
  }

  if (print_slope_area)
  {
    ss.str("");
    ss << name << frame << "_sa";
    slope_area_data( name+"_sa");
  }




}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This routine prints rasters according to some internal model switched
// Added by FJC based on print_rasters function, but prints the metadata as a
// csv which can be read easily by pandas.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::print_rasters_and_csv( int frame )
{
  string outfile_format = "bil";

  cout << endl;
  static ofstream outfile;
  if (outfile.is_open() == false)
  {
    string metadata_fname =  name+"_model_info.csv";
    cout << "Name of raster metadata file is: " <<  metadata_fname << endl;
          outfile.open(metadata_fname.c_str());
          outfile << "Frame_num,";
    outfile << "Time,";
    outfile << "K,";
    outfile << "D,";
    outfile << "Erosion,";
    outfile << "Max_uplift";
    outfile << endl;
  }
  outfile << frame << ",";
  outfile << current_time << ",";
  outfile << get_K() << ",";
  outfile << get_D() << ",";
  outfile << erosion << ",";
  outfile << get_max_uplift() ; // The last comma was generating bug here
  outfile << endl;
  cout << "UPDATE_WARNING::I removed an extra comma from csv file here. It was creating a bug with pandas when reading csv (the python package, not the bamboo junkies). Let me know if it impacts your new outputs" << endl;
  map<string,string> GRS = get_GeoReferencingStrings();

  //cout << "Printing, print elevation is " << print_elevation
  //     << " and erosion is " << print_erosion << endl;

  stringstream ss;
  if (print_elevation)
  {
    if(frame >0 )
      ss << name << frame;
    else
      ss << name << "_init";
    this->write_raster(ss.str(), outfile_format);
  }
  if (print_hillshade)
  {
    cout << "Printing the hillshade" << endl;
    ss.str("");
    ss << name << frame << "_hs";
    LSDRaster * hillshade;
    hillshade = new LSDRaster(*this);
    *hillshade = this->hillshade(45, 315, 1);
    hillshade->write_raster(ss.str(), outfile_format);
    delete hillshade;
  }
  if (print_erosion)
  {
    ss.str("");
    ss << name << frame << "_erosion";
    Array2D <float> erosion_field = calculate_erosion_rates( );
    LSDRaster * erosion;
    erosion = new LSDRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, erosion_field,GRS);
    erosion->write_raster(ss.str(), outfile_format);
    delete erosion;
  }

  if (print_slope_area)
  {
    ss.str("");
    ss << name << frame << "_sa";
    slope_area_data( name+"_sa");
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function cleans the static outfiles
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::close_static_outfiles()
{
  static ofstream er_outfile;
  static ofstream outfile;
  if (er_outfile.is_open())
  {
    er_outfile.close();
  }
  if (outfile.is_open())
  {
    outfile.close();
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates slope-area data and prints to file
// It is retainaed from JAJ's original so that there are no errors with the new
// function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::slope_area_data( string name )
{
  ofstream outfile;
  outfile.open((name).c_str());

  LSDRaster slope;
  Array2D<float> a, b, c, d, e, f;
  Array2D<float> drainage_array(NRows, NCols, 0.0);
  LSDFlowInfo flowData(boundary_conditions, *this);
  int node;
  float DR2 = DataResolution*DataResolution;

  calculate_polyfit_coefficient_matrices(DataResolution, a, b, c, d, e, f);
  slope = calculate_polyfit_slope(d, e);

  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      node = flowData.retrieve_node_from_row_and_column(i, j);
      drainage_array[i][j] = flowData.retrieve_contributing_pixels_of_node(node) * DR2;
    }
  }
  LSDRaster drainage(NRows, NCols, XMinimum, YMinimum, DataResolution,
                     NoDataValue, drainage_array);

  drainage.calculate_polyfit_coefficient_matrices(DataResolution, a, b, c, d, e, f);
  drainage = drainage.calculate_polyfit_elevation(f);

  // Write data
  outfile << name << endl;
  outfile << "Elevation\tSlope\tArea" << endl;

  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      if (slope.get_data_element(i,j) == NoDataValue || RasterData[i][j] == NoDataValue)
      {
        continue;
      }
      else
      {
        outfile << RasterData[i][j];
        outfile << "\t" << slope.get_data_element(i, j);
        outfile << "\t" << drainage.get_data_element(i, j);
        outfile << endl;
      }
        }
  }

  outfile.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print slope area data
// This is an overloaded function that calcualtes slope area
// data based on flags. There are two flag, one for the slope
// calculation and one for the area calculation
// slope_flag is a flag for calculation of the topographic slope
//   0 == polyfit using the data resolution as the smoothing diameter
//   1 == slope calculated with the D8 slopes, with dx = data resolution
//        or DataResolution*sqrt(2) depending on flow direction.
// area_flag is a flag for calculation of the area
//   0 == area using contributing pixels but smoothed to data resolution with polyfit
//   1 == area using contributing pixels only
// The data is printed out in tab seperated format,
// elevation  slope area  predicted_slope
// SMM
// 18/06/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::slope_area_data( string fname, int slope_flag, int area_flag  )
{
  if ( quiet == false )
  {
    cout << "Printing slope-area data, filename is " << fname << endl;
  }

  // open the file for printing
  ofstream outfile;
  outfile.open((fname).c_str());

  // get the flow info
  LSDFlowInfo flowData(boundary_conditions, *this);

  // first calculate the slope
  // slope_flag is a flag for calculation of the topographic slope
  //   0 == polyfit using the data resolution as the smoothing diameter
  //   1 == slope calculated with the D8 slopes, with dx = data resolution
  //        or DataResolution*sqrt(2) depending on flow direction.
  Array2D<float> slope_array(NRows, NCols, 0.0);
  LSDRaster slope;

  if (quiet == false)
  {
    cout << " LSDRasterModel::slope_area_data, slope_flag is: " << slope_flag
         << " and area_flag is: " << area_flag << endl;
  }

  switch (slope_flag)
  {
    case 0:
    {
      vector<LSDRaster> surface_fitting;
      vector<int> raster_selection(8, 0);
      raster_selection[1] = 1;            // this indicates you only want the slope
      surface_fitting = calculate_polyfit_surface_metrics(DataResolution+0.001, raster_selection);
      slope = surface_fitting[1];
    }
    break;
    case 1:
    {
      int n_FI_nodes = flowData.get_NDataNodes();
      int fl_code;
      float donor_elevation;
      float receiver_elevation;
      float dx = DataResolution;
      float dx_rt2 = DataResolution*sqrt(2.0);
      int row, col, rnode, rrow, rcol;
      // loop through the nodes, collecting slope information
      for(int node = 0; node< n_FI_nodes; node++)
      {
         // get the donor elevation
         flowData.retrieve_current_row_and_col(node,row, col);
         donor_elevation =  RasterData[row][col];

         // now the receiver elevation
         flowData.retrieve_receiver_information(node,rnode, rrow, rcol);
         receiver_elevation = RasterData[rrow][rcol];

         // now get the flow length code
         fl_code = flowData.retrieve_flow_length_code_of_node(node);

         // now calculate the slope
         if (fl_code == 1)
         {
           slope_array[row][col] = (donor_elevation-receiver_elevation)/dx;
         }
         else if (fl_code == 2)
         {
           slope_array[row][col] = (donor_elevation-receiver_elevation)/dx_rt2;
         }
      }       // end loop through nodes
      LSDRaster slope2(NRows, NCols, XMinimum, YMinimum, DataResolution,
                         NoDataValue, slope_array);
      slope = slope2;       // a bit of a stupid way to do this but don't want
                            // local slope raster
    }
    break;
    default:
      cout << "Slope_area function, you have chose an incorrect slope flag" << endl;
      exit(EXIT_FAILURE);
  }           // end slope flag switch logic


  // now calculate the area
  // 0 == drainage area using D8 but then smooth the area using polyfit
  // 1 == drainage area D8 No Smoothing
  LSDRaster drainage_area;
  float DR2 = DataResolution*DataResolution;
  Array2D<float> drainage_array(NRows, NCols, 0.0);
  switch (area_flag)
  {
    case 0:
    {
      vector<LSDRaster> surface_fitting_area;
      vector<int> raster_selection(8, 0);
      raster_selection[0] = 1;            // this indicates you only want averaged selection
      int node;

      // loop through the nodes collecting drainage area
      for (int i=0; i<NRows; ++i)
      {
        for (int j=0; j<NCols; ++j)
        {
          node = flowData.retrieve_node_from_row_and_column(i, j);
          drainage_array[i][j] = flowData.retrieve_contributing_pixels_of_node(node) * DR2;
        }
      }
      LSDRaster drainage(NRows, NCols, XMinimum, YMinimum, DataResolution,
                         NoDataValue, drainage_array);

      // fit the surface
      surface_fitting_area = drainage.calculate_polyfit_surface_metrics(DataResolution+0.5*DataResolution, raster_selection);
      drainage_area = surface_fitting_area[0];
    }
    break;
    case 1:
    {
      int node;
      // loop through the nodes getting drainage area.
      for (int i=0; i<NRows; ++i)
      {
        for (int j=0; j<NCols; ++j)
        {
          node = flowData.retrieve_node_from_row_and_column(i, j);
          drainage_array[i][j] = flowData.retrieve_contributing_pixels_of_node(node) * DR2;
        }
      }
      LSDRaster drainage(NRows, NCols, XMinimum, YMinimum, DataResolution,
                         NoDataValue, drainage_array);
      drainage_area = drainage;
    }
    break;
    default:
      cout << "Slope_area function, you have chose an incorrect area flag" << endl;
       exit(EXIT_FAILURE);
  }       // end area flag switch logic

  // Write data
  outfile << name << endl;
  outfile << "row\tcol\tElevation\tSlope\tArea\tPredicted_slope" << endl;

  float predicted_slope, A_pow, U_over_KApowm;
  int n_sa_nodes = 0;
  float err;
  float err_tot = 0;

  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      if (slope.get_data_element(i,j) == NoDataValue
            || RasterData[i][j] == NoDataValue
            || is_base_level(i,j))
      {
        continue;
      }
      else
      {
        A_pow = pow( drainage_area.get_data_element(i, j), m);

        if (A_pow == 0)
        {
          predicted_slope = NoDataValue;
        }
        else
        {
          U_over_KApowm = get_uplift_rate_at_cell(i,j)/(get_K()*A_pow);
          predicted_slope = pow(U_over_KApowm, (1/n));
        }

        err = (slope.get_data_element(i, j)- predicted_slope)
              *(slope.get_data_element(i, j)- predicted_slope)
              /  (predicted_slope*predicted_slope);
        err_tot+= err;
        n_sa_nodes++;

        outfile << i << "\t" << j << "\t" << RasterData[i][j];
        outfile << "\t" << slope.get_data_element(i, j);
        outfile << "\t" << drainage_area.get_data_element(i, j);
        outfile << "\t" << predicted_slope;
        outfile << endl;
      }
    }
  }

  if (quiet == false)
  {
    cout << "mean error % between predicted and measured slope is: " << err_tot/n_sa_nodes << endl;
  }
  outfile.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes the default parameters and prints out a template
// paramter file. It is used so that if you run the model with
// no arguments you can look back at what parameters were used.
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::make_template_param_file(string filename)
{
  ofstream param;
  param.open(filename.c_str());

  param << "# Template for parameter file" << endl;
  param << "Run Name:\t\ttemplate" << endl;
  param << "NRows:\t\t\t100" << endl;
  param << "NCols:\t\t\t100" << endl;
  param << "Resolution:\t\t1" << endl;
  param << "Boundary code:\t\tbnbn\tNorth, east, south, west" << endl;
  param << "# b = base level, p = periodic, n = no flow (default)" << endl;
  param << "Time step:\t\t50" << endl;
  param << "End time:\t\t2000" << endl;
  param << "End time mode:\t\t0\t(if 1, wait for steady state to set the time to count down)" << endl;
  param << "Uplift mode:\t\t0\tBlock uplift" << endl;
  param << "Max uplift:\t\t0.001" << endl;
  param << "Tolerance:\t\t0.0001" << endl;
  param << "Print interval:\t\t5" << endl;
  param << "#Periodicity:\t\t1000" << endl;

  param << "\n#####################" << endl;
  param << "Fluvial:\t\ton" << endl;
  param << "K:\t\t\t0.01" << endl;
  param << "m:\t\t\t0.5" << endl;
  param << "n:\t\t\t1" << endl;
  param << "K mode:\t\t\t0\tconstant" << endl;
  param << "#K amplitude:\t\t0.005" << endl;

  param << "\n#####################" << endl;
  param << "Hillslope:\t\ton" << endl;
  param << "Non-linear:\t\toff" << endl;
  param << "Threshold drainage:\t-1\t(if negative, ignored)" << endl;
  param << "D:\t\t\t0.05" << endl;
  param << "S_c:\t\t\t30\tdegrees" << endl;
  param << "D mode:\t\t\t0\tConstant" << endl;
  param << "#D amplitude:\t\t0.005" << endl;

  param << "\n#####################" << endl;
  param << "Isostasy:\t\toff" << endl;
  param << "Flexure:\t\toff" << endl;
  param << "Rigidity:\t\t1000000" << endl;

  param.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the fluvial erodability. It has a number of switches
// that determine how K is calcualted.
// K_mode == 1 sine wave
// K_mode == 2 square wave
// K_mode == 3 read from file
// K_mode: default is constant value
//
// Author JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_K( void )
{

  static bool copied = false;

  if (K_mode == 3 && copied == false)
  {
    stringstream ss;
    ss << ".K_file_" << name << ".aux";
    ifstream infile("K_file", ios::binary);
    ofstream outfile(ss.str().c_str(), ios::binary);

    outfile << infile.rdbuf();

    // This inhibits portablity
    // SMM: Something JAJ put in to manage permissions but I've got rid of it.
    //ss.str("");
    //ss << "chmod -w .K_file_" << name << ".aux";
    //system(ss.str().c_str());

    copied = true;
  }

  // in all of these switches, the ternary operator is used: a constant value
  // is selected if initial_steady_state is false, otherwise the parameter is
  // determined from a function
  switch( K_mode ) {
    case 1:      // sin wave
      return (initial_steady_state || cycle_steady_check) ? periodic_parameter( K_fluv, K_amplitude ) : K_fluv;
      break;
    case 2:      // square wave
      return (initial_steady_state) ? square_wave_parameter( K_fluv, K_amplitude ) : K_fluv;
      break;
    case 3:
      return (initial_steady_state) ? stream_K_fluv() : K_fluv;
      break;
    default:    // constant
      return K_fluv;
      break;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the soil transport coefficient. It has a number of switches
// that determine how D is calcualted.
// D_mode == 1 sine wave
// D_mode == 2 square wave
// D_mode == 3 read from file
// D_mode: default is constant value
//
// Author JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::get_D( void )
{
  // in all of these switches, the ternary operator is used: a constant value
  // is selected if initial_steady_state is false, otherwise the parameter is
  // determined from a function

  //cout << "LINE 454 Getting D, D mode is: " << D_mode << " and D amp is: " << D_amplitude << endl;
  float thisD;
  //cout << "initial steady state: " << initial_steady_state << endl;
  static bool copied = false;

  if (D_mode == 3 && copied == false)
  {
    stringstream ss;
    ss << ".D_file_" << name << ".aux";
    ifstream infile("D_file", ios::binary);
    ofstream outfile(ss.str().c_str(), ios::binary);

    outfile << infile.rdbuf();

    // This inhibits portablity
    // SMM: Something JAJ put in to manage permissions but I've got rid of it.
    //ss.str("");
    //ss << "chmod -w .D_file_" << name << ".aux";
    //system(ss.str().c_str());

    copied = true;
  }

  switch( D_mode )
  {
    case 1:
    {
      // sin wave
      thisD =  ((initial_steady_state) ? periodic_parameter( K_soil, D_amplitude ) : K_soil);
      //cout << "LINE 5568 Sine wave: D is: " << thisD << " K_soil is: " << K_soil << endl;
      break;
    }
    case 2:   // square wave
    {
      thisD =  ((initial_steady_state) ? square_wave_parameter( K_soil, D_amplitude ) : K_soil);
      break;
    }
    case 3:
    {
      thisD = ((initial_steady_state) ? stream_K_soil() : K_soil);
      break;
    }
    default:    // constant
    {
      thisD =  K_soil;
      break;
    }
  }
  //cout << "thisD is: " << thisD << endl;
  return thisD;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates the value of a sinusoidal periodic variable.
// To calcualte the variable, it uses the data member current_time
// and delay_switch to calcualte where in the cycle the parameter is.
// It uses a 'period mode' variable to determine the nature of the periodic
// forcing (which I'll have to look up later, SMM)
// Author JAJ sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::periodic_parameter( float base_param, float amplitude )
{
  float result;

  if (period_mode == 3 || period_mode == 4)
  {
    result = p_weight * sin( (current_time-time_delay-switch_delay)*2*PI/periodicity )*amplitude +
        (1-p_weight) * sin( (current_time-time_delay-switch_delay)*2*PI/periodicity_2)*amplitude + base_param;
  }
  else
  {
    //cout << "LINE 5611: ctime: " << current_time << " td: " << time_delay << " switch_delay: " << switch_delay << endl;
    //cout << "amplitude: " << amplitude << " base_param: " << base_param << " periodicity: " << periodicity << endl;
    result = sin( (current_time - time_delay - switch_delay) * 2 * PI / periodicity )* amplitude + base_param;
  }

  return result;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates the value of a sinusoidal periodic variable.
// To calcualte the variable, it uses the data member current_time
// and delay_switch to calcualte where in the cycle the parameter is.
// Author JAJ sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::square_wave_parameter( float base_param, float amplitude )
{
  int wave = (current_time - time_delay - switch_delay) /  (this->periodicity / 2);
  if (wave % 2 == 0)   wave =  1;
  else                 wave = -1;

  return base_param + (wave*amplitude);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function reads in the K parameter as a stream
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::stream_K_fluv( void )
{
  static float upr_param = K_fluv;
  static float lwr_param = K_fluv;
  static float upr_t = -99;
  static float lwr_t = 0;
  static ifstream strm;
  if (strm.is_open() == false)
  {
    stringstream ss;
    ss << ".K_file_" << name << ".aux";
    strm.open(ss.str().c_str());
  }

  float temp;
  bool read = true;

  while (current_time >= upr_t)
  {
    if (strm >> temp)
    {
      if (upr_t == -99)
        lwr_t = time_delay;
      else
        lwr_t = upr_t;
      lwr_param = upr_param;
      upr_t = temp+time_delay;
      strm >> upr_param;
      read = true;
    }
    else
    {
      read = false;
      break;
    }

  }
  if (read)
    return (upr_param - lwr_param) * (current_time-lwr_t) / (upr_t - lwr_t) + lwr_param;
  else
    return upr_param;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function reads in the D parameter as a stream
// JAJ, sometime 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDRasterModel::stream_K_soil( void )
{
  static float upr_param = K_soil;
  static float lwr_param = K_soil;
  static float upr_t = -99;
  static float lwr_t = 0;
  static ifstream strm;
  if (strm.is_open() == false)
  {
    stringstream ss;
    ss << ".D_file_" << name << ".aux";
    strm.open(ss.str().c_str());
  }

  float temp;
  bool read;

  while (current_time >= upr_t)
  {
    if (strm >> temp)
    {
      if (upr_t == -99)
        lwr_t = time_delay;
      else
        lwr_t = upr_t;
      lwr_param = upr_param;
      upr_t = temp+time_delay;
      strm >> upr_param;
      read = true;
    }
    else
      read = false;

  }
  if (read)
    return (upr_param - lwr_param) * (current_time-lwr_t) / (upr_t - lwr_t) + lwr_param;
  else
    return upr_param;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This snaps the period to a timestep
// its main use is to make sure the max and min paramter values are reached
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::snap_periodicity( void )
{
  periodicity = ceil(periodicity/timeStep) * timeStep;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets data members associated with the MuddPILE solution of
// the nonlinear hillslope equations
// SMM, 01/07/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::MuddPILE_initiate_assembler_matrix(void)
{
  float dx = DataResolution;
  float D_nl = get_D();

  //cout << "LSDRasterModel::MuddPILE_initiate_assembler_matrix, D_nl is: "
  //     << D_nl << " and S_c is: " << S_c << endl;

  // update data_members
  inv_dx_S_c_squared = 1/(dx*dx*S_c*S_c);
  dx_front_term = timeStep*D_nl/(dx*dx);
  problem_dimension = NRows*NCols;

  //cout << "MuddPILE_initiate_assembler_matrix, calling k values" << endl;
  // this sets the vector k data members
  MuddPILE_calculate_k_values_for_assembly_matrix();

  //cout << "LSDRasterModel::MuddPILE_initiate_assembler_matrix(void) initiated! " << endl;

  // check that zeta_last_iter has been initiated
  //if(zeta_last_iter.dim1() != NRows || zeta_last_iter.dim2() != NCols)
  //{
  //  cout << "LSDRasterModel::MuddPILE_initiate_assembler_matrix, Warning " << endl
  //       << "zeta_last_iter doesn't exist" << endl;
  // }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This calculates the k values needed for the assembly matrix
// this function creates vectors of integers that refer to the k values, that is
// the index into the vectorised matrix of zeta values, that is used in the assembly matrix
// the number of elements in the k vectors is N_rows*N_cols
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::MuddPILE_calculate_k_values_for_assembly_matrix()
{
  // SMM: need to double check if this is only true of periodic boundary
  // conditions. It seems that it is the case.
  int N_elements_in_k_vec = NRows*NCols;

  // initialise the vectors with empty values
  vector<int> empty_vec(N_elements_in_k_vec,0);
  vec_k_value_i_j   = empty_vec;
  vec_k_value_ip1_j = empty_vec;
  vec_k_value_im1_j = empty_vec;
  vec_k_value_i_jp1 = empty_vec;
  vec_k_value_i_jm1 = empty_vec;

  //cout << "MuddPILE calculate k vecs, n elements: " << N_elements_in_k_vec << endl;

  // we loop through each node
  // This is done without buffering.
  // Some of these entries will result in negative indices, but these
  // will never be called since the boundaries are set by the
  // constant elevation boundary condition
  int counter = 0;
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      // bounds checking
      if (counter > N_elements_in_k_vec)
      {
         cout << "DANGER, LSDRasterModel::MuddPILE_calculate_k_values_for_assembly_matrix" << endl
              << "counter is out of bounds!" << endl;
       }


      vec_k_value_ip1_j[counter] = NCols*(row+1)+col;
      vec_k_value_im1_j[counter] = NCols*(row-1)+col;
      vec_k_value_i_j[counter] = NCols*(row)+col;

      // logic for west periodic boundary
      if(col == 0)
      {
        vec_k_value_i_jp1[counter] = NCols*(row)+col+1;
        vec_k_value_i_jm1[counter] = NCols*(row)+NCols-1;
      }
      // logic for east periodic boundary
      else if(col == NCols-1)
      {
        vec_k_value_i_jp1[counter] = NCols*(row);
        vec_k_value_i_jm1[counter] = NCols*(row)+col-1;
      }
      // logic for rest of matrix
      else
      {
        vec_k_value_i_jp1[counter] = NCols*(row)+col+1;
        vec_k_value_i_jm1[counter] = NCols*(row)+col-1;
      }

      // increment counter
      counter++;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function assembles the solution matrix for nonlinear creep transport.
// There is no buffereing of the surface: the constant elevation nodes are
// in the zeta value within the domain of the DataArray
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::MuddPILE_assemble_matrix(Array2D<float>& uplift_rate,
             Array2D<float>& fluvial_erosion_rate,
             mtl::compressed2D<float>& mtl_Assembly_matrix,
             mtl::dense_vector<float>& mtl_b_vector)
{

  // get the soil diffusivity of the current step
  // note that in the future we might want
  // i) spatially varying D values
  // ii) time and space varying S_c
  float D_nl = get_D();

  //cout<< "LINE 5494 D is: " << D_nl << " zti: " << zeta_this_iter[10][10]
  //    << " zeta_lts: "<< zeta_last_timestep[10][10] << " rd: " << RasterData[10][10] << endl;

  if(D_mode == 1)
  {
    //cout << "Variable D!, D is: " << D_nl << endl;
  }
  dx_front_term = timeStep*D_nl/(DataResolution*DataResolution);

  // get the number of elements in the k vec
  int n_k_elements = vec_k_value_i_j.size();

  // the coefficients in the assembly matrix
  float A,B,C,D;

  // reset the assembly and b vector
  mtl_Assembly_matrix = 0.0;
  mtl_b_vector = 0.0;

  // create the inserter. This is deleted when this function is exited
  mtl::mat::inserter< mtl::compressed2D<float> > ins(mtl_Assembly_matrix);

  // first we assemble the boundary nodes. First the nodes in row 0
  //cout << "Line 5180, getting south boundary" << endl;
  for (int k = 0; k<NCols; k++)
  {
    //cout << "k is " << k << endl;

    ins[k][k] << 1.0;
    //cout << "inserted boundary, now doing b vec" << endl;
    //cout << "b vec this k" << mtl_b_vector[k] << endl;
    //cout << "zeta_last_iter[0][0] " << RasterData[0][0] << endl;
    mtl_b_vector[k] =  RasterData[0][0];
    //cout << "did b vec" << endl;
  }


  // now assemble the north boundary
  // in this implementation there is no buffered surface
  int starting_north_boundary = (NRows-1)*(NCols);
  int one_past_last_north_boundary = (NRows)*NCols;
  //cout << "Line 5190, getting N boundary!" << endl;
  for (int k = starting_north_boundary; k < one_past_last_north_boundary; k++)
  {
    ins[k][k] << 1.0;
    mtl_b_vector[k] = RasterData[NRows-1][0];
  }



  // now assemble the rest
  // we loop through each node
  int counter = NCols;       // the counter starts at NCols because the assumbly
                             // matrix starts at the first row.
  float b_value;
  int k_value_i_j,k_value_ip1_j,k_value_im1_j,k_value_i_jp1,k_value_i_jm1;

  // this does not loop over row 0 or NRows-1 because these have a
  // fixed elevation, these are addressed earlier in the boundary nodes
  //cout << "Line 5204, assembling!!!" << endl;
  for (int row = 1; row<NRows-1; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      // bounds check
      if (counter >= n_k_elements)
      {
        cout << "Danger!!!, counter: " << counter << " n_k: " << n_k_elements << endl;
      }

      b_value = zeta_last_timestep[row][col]+timeStep*uplift_rate[row][col]
                                   -timeStep*fluvial_erosion_rate[row][col];

      k_value_ip1_j = vec_k_value_ip1_j[counter];
      k_value_im1_j = vec_k_value_im1_j[counter];
      k_value_i_j   = vec_k_value_i_j[counter];
      k_value_i_jp1 = vec_k_value_i_jp1[counter];
      k_value_i_jm1 = vec_k_value_i_jm1[counter];

      // bounds check
     if (k_value_ip1_j >= problem_dimension)
      {
        cout << "Warning, k i+1,j value out of bounds!" << endl;
        cout << "problem dimension: " << problem_dimension << endl;
        cout << "value: " <<  k_value_ip1_j << endl;
        cout << "row" << row << " and col: " << col << endl;
      }
      if (k_value_im1_j >= problem_dimension)
      {
        cout << "Warning, k i-1,j value out of bounds!" << endl;
        cout << "problem dimension: " << problem_dimension << endl;
        cout << "value: " <<  k_value_im1_j << endl;
        cout << "row" << row << " and col: " << col << endl;
      }
      if (k_value_i_j >= problem_dimension)
      {
        cout << "Warning, k i,j value out of bounds!" << endl;
        cout << "problem dimension: " << problem_dimension << endl;
        cout << "value: " <<  k_value_i_j << endl;
        cout << "row" << row << " and col: " << col << endl;
      }
      if (k_value_i_jm1 >= problem_dimension)
      {
        cout << "Warning, k i,j-1 value out of bounds!" << endl;
        cout << "problem dimension: " << problem_dimension << endl;
        cout << "value: " <<  k_value_i_jm1 << endl;
        cout << "row" << row << " and col: " << col << endl;
      }
      if (k_value_i_jp1 >= problem_dimension)
      {
        cout << "Warning, k i,j+1 value out of bounds!" << endl;
        cout << "problem dimension: " << problem_dimension << endl;
        cout << "value: " <<  k_value_i_jp1 << endl;
        cout << "row" << row << " and col: " << col << endl;
      }

      A =  dx_front_term/(1 -
              (RasterData[row+1][col]-RasterData[row][col])*
              (RasterData[row+1][col]-RasterData[row][col])*
          inv_dx_S_c_squared);
      B = dx_front_term/(1 -
              (RasterData[row][col]-RasterData[row-1][col])*
              (RasterData[row][col]-RasterData[row-1][col])*
          inv_dx_S_c_squared);

      // logic for west periodic boundary
      if(col == 0)
      {
        C = dx_front_term/(1 -
                 (RasterData[row][col+1]-RasterData[row][col])*
                 (RasterData[row][col+1]-RasterData[row][col])*
             inv_dx_S_c_squared);
        D = dx_front_term/(1 -
                 (RasterData[row][col]-RasterData[row][NCols-1])*
                 (RasterData[row][col]-RasterData[row][NCols-1])*
             inv_dx_S_c_squared);
      }
      // logic for east periodic boundary
      else if(col == NCols-1)
      {
        C = dx_front_term/(1 -
                 (RasterData[row][0]-RasterData[row][col])*
                 (RasterData[row][0]-RasterData[row][col])*
             inv_dx_S_c_squared);

        D = dx_front_term/(1 -
                 (RasterData[row][col]-RasterData[row][col-1])*
                 (RasterData[row][col]-RasterData[row][col-1])*
             inv_dx_S_c_squared);

      }
      // logic for rest of matrix
      else
      {
        C = dx_front_term/(1 -
                 (RasterData[row][col+1]-RasterData[row][col])*
                 (RasterData[row][col+1]-RasterData[row][col])*
             inv_dx_S_c_squared);

        D = dx_front_term/(1 -
                 (RasterData[row][col]-RasterData[row][col-1])*
                 (RasterData[row][col]-RasterData[row][col-1])*
             inv_dx_S_c_squared);

      }


      // some couts for bug checking
      //cout << "problem dimension: " << problem_dimension << ", k values" << endl;
      //cout << "row: " << row << " and col " << col << endl;
      //cout << "k i,j: " <<  k_value_i_j << endl;
      //cout << "k i+1,j: " <<  k_value_ip1_j << endl;
      //cout << "k i-1,j: " <<  k_value_im1_j << endl;
      //cout << "k i,j+1: " <<  k_value_i_jp1 << endl;
      //cout << "k i,j-1: " <<  k_value_i_jm1 << endl;

      // place the values in the assembly matrix and the b vector
      mtl_b_vector[k_value_i_j] = b_value;
      ins[k_value_i_j][k_value_ip1_j] << -A;
      ins[k_value_i_j][k_value_im1_j] << -B;
      ins[k_value_i_j][k_value_i_jp1] << -C;
      ins[k_value_i_j][k_value_i_jm1] << -D;
      ins[k_value_i_j][k_value_i_j] << 1+A+B+C+D;

      counter++;
    }
  }

  //cout << "Line 6580 assembled matrix " << endl;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function assembles the solution matrix
// It used the mtl library for sparse matrices
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::MuddPILE_solve_assembler_matrix(Array2D<float>& uplift_rate,
                        Array2D<float>& fluvial_erosion_rate)
{

  // reset the zeta array for this iteration
  Array2D<float> empty_zeta(NRows,NCols,0.0);
  zeta_this_iter = empty_zeta.copy();

  // create a mtl matrix
  // NOTE: you could probably save time by creating the mtl matrix and vector
  // in main()
  mtl::compressed2D<float> mtl_Assembly_matrix(problem_dimension, problem_dimension);
  mtl::dense_vector<float> mtl_b_vector(problem_dimension,0.0);

  //cout << "dense vector is: " << mtl_b_vector << endl;

  //cout<< "LINE 5701 zti: " << zeta_this_iter[10][10]
  //    << " zeta_lts: " << zeta_last_timestep[10][10] << " rd: " << RasterData[10][10] << endl;

  /*
  if(current_time >= 11000)
  {
     cout <<"YO 5806\n";
      string asc_name = "asc";
      string this_time = itoa(int(current_time));
      string RDfname = "RDassem_t"+this_time;
      LSDRaster RD(NRows, NCols, XMinimum, YMinimum,
            DataResolution, NoDataValue, RasterData.copy());
      RD.write_raster(RDfname,asc_name);

      string ZTIfname = "ZTIassem_t"+this_time;
      LSDRaster ZTI(NRows, NCols, XMinimum, YMinimum,
            DataResolution, NoDataValue, zeta_this_iter.copy());
      ZTI.write_raster(ZTIfname,asc_name);

      string LTSfname = "LTSassem_t"+this_time;
      LSDRaster LTS(NRows, NCols, XMinimum, YMinimum,
            DataResolution, NoDataValue, RasterData.copy());
      LTS.write_raster(LTSfname,asc_name);
    }
  */

  // assemble the matrix
  //cout << "LINE 5289 assembling matrix 1st time" << endl;
  //cout << "Line 5309, problem dimension: " << problem_dimension << endl;
  MuddPILE_assemble_matrix(uplift_rate, fluvial_erosion_rate,mtl_Assembly_matrix,
                           mtl_b_vector);
  //cout << "LINE 5292 assembled!" << endl;


  //cout<< "LINE 5714 zti: " << zeta_this_iter[10][10]
  //    << " zeta_lts: " << zeta_last_timestep[10][10] << " rd: " << RasterData[10][10] << endl;
  //cout<< "LINE 5716 zti: [0][10] " << zeta_this_iter[0][10]
  //    << " zeta_lts: "<< zeta_last_timestep[0][10] << " rd: " << RasterData[0][10] << endl;

  // some couts for bug checking
  //cout << "matrix assembled!" << endl;
  //ofstream assembly_out;
  //assembly_out.open("assembly.data");
  //assembly_out << mtl_Assembly_matrix << endl;
  //assembly_out.close();

  // now solve the mtl system
  // Create an ILU(0) preconditioner
  long time_start, time_end, time_diff;
  bool show_time = false;
  time_start = time(NULL);
  //cout << "LINE 5404 making P " << endl;
  itl::pc::ilu_0< mtl::compressed2D<float> > P(mtl_Assembly_matrix);
  //cout << "LINE 5405 made P " << endl;
  mtl::dense_vector<float> mtl_zeta_solved_vector(problem_dimension);
  //cout << "LINE 5406 made mtl_zeta_solved_vector" << endl;
  itl::basic_iteration<float> iter(mtl_b_vector, 500, 1.e-8);
  //cout << "LINE 5408 made iter" << endl;
  bicgstab(mtl_Assembly_matrix, mtl_zeta_solved_vector, mtl_b_vector, P, iter);
  time_end = time(NULL);
  time_diff = time_end-time_start;

  if(show_time)
  {
    cout << "iter MTL bicg took: " << time_diff << endl;
  }

  // now reconstitute zeta
  int counter = 0;
  //cout << "LINE 5413 reconstituting data!" << endl;
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      //cout << "counter is: " << counter << endl;
      zeta_this_iter[row][col] = mtl_zeta_solved_vector[counter];
      counter++;
    }
  }

  /*
  if(current_time >= 11000)
    {
      cout <<"YO 5806\n";
      string asc_name = "asc";
      string this_time = itoa(int(current_time));
      string RDfname = "RDinterm_t"+this_time;
      LSDRaster RD(NRows, NCols, XMinimum, YMinimum,
            DataResolution, NoDataValue, RasterData.copy());
      RD.write_raster(RDfname,asc_name);

      string ZTIfname = "ZTIinterm_t"+this_time;
      LSDRaster ZTI(NRows, NCols, XMinimum, YMinimum,
            DataResolution, NoDataValue, zeta_this_iter.copy());
      ZTI.write_raster(ZTIfname,asc_name);

      string LTSfname = "LTSinterm_t"+this_time;
      LSDRaster LTS(NRows, NCols, XMinimum, YMinimum,
            DataResolution, NoDataValue, RasterData.copy());
      LTS.write_raster(LTSfname,asc_name);
    }
   */
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// do a creep timestep
// NOTE you need to run MuddPILE_initiate_assembler_matrix before you run this function
// At the end of this iteration RasterData will have the new surface elevations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::MuddPILE_nonlinear_creep_timestep(Array2D<float>& uplift_rate,
            Array2D<float>& fluvial_erosion_rate,
            float iteration_tolerance)
{

  // check to see if the model has been initiated
  if ( int(vec_k_value_i_j.size()) != (NRows*NCols) )
  {
    cout << "LSDRasterModel::MuddPILE creep, initialising k vectors for 1st time" << endl;
    MuddPILE_initiate_assembler_matrix();
  }

  // make sure that the N and S boundaries are at zero
  for (int col = 0; col<NCols; col++)
  {
    RasterData[0][col] = 0;
    RasterData[NRows-1][col] = 0;

    //zeta_last_timestep[0][col] = 0;
    //zeta_last_timestep[NRows-1][col] = 0;

  }

  // reset old zeta
  zeta_last_timestep = RasterData.copy();

  // reset the zeta_this_iter
  zeta_this_iter = RasterData.copy();

  // set up residual
  float residual;
  float N_nodes = float(NRows*NCols);
  int iteration = 0;
  int Max_iter = 100;
  do
  {
    residual = 0.0;

    //cout << "Time is: " << current_time << endl;
    //cout << "LINE 5775 zti[10][10]: " << zeta_this_iter[10][10] << " and uplift: "
    //     << uplift_rate[10][10] << " and fluv: " << fluvial_erosion_rate[10][10] << endl;

    // this solves for zeta_this_iter
    MuddPILE_solve_assembler_matrix(uplift_rate, fluvial_erosion_rate);

    //cout << "LINE 5784 zti[10][10]: " << zeta_this_iter[10][10] << endl;

    // check the residuals (basically this is the aveage elevation change between intermediate
    // zeta values
    for (int row = 0; row<NRows; row++)
    {
      for (int col = 0; col<NCols; col++)
      {
        // in the first iteration, RasterData contains the elevations from the
        // previous timestep. In subsequent iterations it contains the last iteration
        residual+= sqrt( (zeta_this_iter[row][col]-RasterData[row][col])*
                 (zeta_this_iter[row][col]-RasterData[row][col]) );
      }
    }
    residual = residual/N_nodes;

    // reset the last iteration of surface elevations zeta
    RasterData = zeta_this_iter.copy();
    iteration++;

    if (iteration%5 == 0)
    {
      std::cout << "iteration is: " << iteration << " and residual RMSE is: " << residual << endl;
    }
    if (iteration > Max_iter)
    {
      iteration_tolerance = iteration_tolerance*10;
      iteration = 0;
    }

  } while (residual > iteration_tolerance);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Run the soil diffusion routine, but with a set tolerance
// and no fluvial or tectonic uplift. The end result is an updated RasterData
// surface.
//
// This particular version is set to run within JAJ's 'run_components'
// module
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterModel::MuddPILE_nl_soil_diffusion_nouplift()
{

  // set the fluvial and uplift rasters to zero
  // I don't use the same raster for both since I'm a bit worried about
  // passing the same reference to a data oject for two different computations
  // in the nonlinear_creep_timestep module
  Array2D<float> zero_uplift(NRows,NCols,0.0);
  Array2D<float> zero_fluvial(NRows,NCols,0.0);

  float default_tolerance = 3e-6;

  // run a timestep

  //cout << "Line 5841, data[10][10]: " << RasterData[10][10] << endl;
  MuddPILE_nonlinear_creep_timestep(zero_uplift, zero_fluvial,default_tolerance);
  //cout << "Line 5843, data[10][10]: " << RasterData[10][10] << endl;
}





#endif
