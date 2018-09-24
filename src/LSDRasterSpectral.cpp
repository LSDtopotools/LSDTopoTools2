//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterSpectral
// Land Surface Dynamics StatsTools
//
// An object for manipulating rasters developed for the University of Edinburgh
//  Land Surface Dynamics group topographic toolbox. This is a derivative class
// from LSDRaster, for use specifically with spectral analysis.
//
// These tools have been seperated from the LSDRaster class mainly because
//  they require the FFTW library and are therefore less portable than
//  the standard LSDRaster object.
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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterSpectral.cpp
// cpp file for the LSDRasterSpectral object
// LSD stands for Land Surface Dynamics
// This object perform spectral analysis and is seperate from LSDRaster
// simply because it requires the FFTW package so this can be removed
// from compilation to retain portability
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
// Stuart Grieve, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1		02/04/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <string.h>
#include "TNT/tnt.h"
#include "TNT/jama_lu.h"
#include "TNT/jama_eig.h"
#include "LSDRaster.hpp"
#include "LSDRasterSpectral.hpp"
#include "LSDStatsTools.hpp"
#include "LSDIndexRaster.hpp"
#include "fftw-3.3.4/api/fftw3.h"
using namespace std;
using namespace TNT;
using namespace JAMA;

#ifndef LSDRasterSpectral_CPP
#define LSDRasterSpectral_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Assignment operator
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRasterSpectral& LSDRasterSpectral::operator=(const LSDRasterSpectral& rhs)
{
  if (&rhs != this)
  {
    create(rhs.get_NRows(),rhs.get_NCols(),rhs.get_XMinimum(),rhs.get_YMinimum(),
           rhs.get_DataResolution(),rhs.get_NoDataValue(),rhs.get_RasterData(),
           rhs.get_GeoReferencingStrings());
  }
  return *this;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// the create function. This is default and throws an error
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterSpectral::create()
{
  cout << "LSDRasterSpectral line 63 You need to initialize with a filename!" << endl;
  //exit(EXIT_FAILURE);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// create function that creates a LSDSpectralRaster with
// that is square, and that has dimensions NRow, NCols = 2^raster_order
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterSpectral::create(int raster_order, float cellsize, float ndv)
{
  // get the raster size
  int raster_size = pow(2,raster_order);
  cout << "You are making an LSDSpectralRaster with an size of: " << raster_size << endl;

  // now set the raster properties
  NRows = raster_size;
  NCols = raster_size;
  XMinimum = 0.0;
  YMinimum = 0.0;
  DataResolution = cellsize;
  NoDataValue = ndv;
  Ly = int(pow(2,ceil(log(NRows)/log(2))));
  Lx = int(pow(2,ceil(log(NCols)/log(2))));
  dfx = 1/(DataResolution*float(Lx));
  dfy = 1/(DataResolution*float(Ly));
  NyquistFreq = 1/(2*DataResolution);
  WSS=float(NRows*NCols);
  cout << "Created square raster with NRows = " << NRows << " and Ly: " << Ly << endl;
  cout << "NCols: " << NCols << " and Lx: " << Lx << endl;

  Array2D<float> test_data(NRows,NCols,ndv);

  RasterData = test_data.copy();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this creates a raster using an infile
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterSpectral::create(string filename, string extension)
{
  read_raster(filename,extension);
  Ly = int(pow(2,ceil(log(NRows)/log(2))));
  Lx = int(pow(2,ceil(log(NCols)/log(2)))); 
  dfx = 1/(DataResolution*float(Lx));
  dfy = 1/(DataResolution*float(Ly));      
  NyquistFreq = 1/(2*DataResolution);
  WSS=float(NRows*NCols);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this creates a raster filled with some data values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterSpectral::create(int nrows, int ncols, float xmin, float ymin,
            float cellsize, float ndv, Array2D<float> data)
{
  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;
  Ly = int(pow(2,ceil(log(NRows)/log(2))));
  Lx = int(pow(2,ceil(log(NCols)/log(2)))); 
  dfx = 1/(DataResolution*float(Lx));
  dfy = 1/(DataResolution*float(Ly));   
  NyquistFreq = 1/(2*DataResolution);
  WSS=float(NRows*NCols);
  RasterData = data.copy();

  if (RasterData.dim1() != NRows)
  {
    cout << "dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
  if (RasterData.dim2() != NCols)
  {
    cout << "dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Like the above function, but copies the GeoReferencing
// SMM 2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterSpectral::create(int nrows, int ncols, float xmin, float ymin,
           float cellsize, float ndv, Array2D<float> data, map<string,string> temp_GRS)
{
  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;
  GeoReferencingStrings = temp_GRS;
  Ly = int(pow(2,ceil(log(NRows)/log(2))));
  Lx = int(pow(2,ceil(log(NCols)/log(2)))); 
  dfx = 1/(DataResolution*float(Lx));
  dfy = 1/(DataResolution*float(Ly));   
  NyquistFreq = 1/(2*DataResolution);
  WSS=float(NRows*NCols);
  RasterData = data.copy();

  if (RasterData.dim1() != NRows)
  {
    cout << "LSDRaster line 89 dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
  if (RasterData.dim2() != NCols)
  {
    cout << "LSDRaster line 94 dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function creates an LSDSpectralRaster from an LSDRaster
void LSDRasterSpectral::create(LSDRaster& An_LSDRaster)
{
  NRows = An_LSDRaster.get_NRows();
  NCols = An_LSDRaster.get_NCols();
  XMinimum = An_LSDRaster.get_XMinimum();
  YMinimum = An_LSDRaster.get_YMinimum();
  DataResolution = An_LSDRaster.get_DataResolution();
  NoDataValue = An_LSDRaster.get_NoDataValue();
  GeoReferencingStrings = An_LSDRaster.get_GeoReferencingStrings();
  Ly = int(pow(2,ceil(log(NRows)/log(2))));
  Lx = int(pow(2,ceil(log(NCols)/log(2)))); 
  dfx = 1/(DataResolution*float(Lx));
  dfy = 1/(DataResolution*float(Ly));      
  NyquistFreq = 1/(2*DataResolution);
  WSS=float(NRows*NCols);
  RasterData = An_LSDRaster.get_RasterData();
  
  cout << "Loaded a spectralreaster from another raster. Here are the vitalstatistix." << endl;
  cout << "NRows = " << NRows << endl;
  cout << "NCols = " << NCols  << endl;
  cout << "XMinimum = " << XMinimum  << endl;
  cout << "YMinimum = " << YMinimum  << endl;
  cout << "DataResolution =" << DataResolution  << endl;
  cout << "NoDataValue = " <<  NoDataValue << endl;
  cout << "Ly = " << Ly  << endl;
  cout << "Lx = " << Lx  << endl;
  cout << "dfx = " << dfx   << endl;
  cout << "dfy = " << dfy   << endl;
  cout << "WSS= " << WSS  << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Fourier Helper functions
// These function don't really use the data stored within the object
// but are useful in Fourier analysis so are located here as member functions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function returns the frequencies in the direction of the rows
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDRasterSpectral::get_row_direction_frequencies_unshifted()
{
  vector<float> freq_values(NRows);

  // now get the frequencies based on the cellsize.
  int lastrow = 0;
  for(int row = 0; row<=(NRows-1)/2; row++)
  {
    freq_values[row] = float(row)*DataResolution/(float(NRows));
    //cout << "Row: " << row << endl;
    lastrow = row;
  }
  for(int row = (NRows)/2; row > 0; row--)
  {
    lastrow++;
    freq_values[lastrow] = -float(row)*DataResolution/(float(NRows));
    //cout << "Row is: " << lastrow << "and selection is: " << -row << endl;
  }

  return freq_values;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function returns the frequencies in the direction of the cols
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDRasterSpectral::get_col_direction_frequencies_unshifted()
{
  vector<float> freq_values(NCols);

  // now get the frequencies based on the cellsize.
  int lastcol = 0;
  for(int col = 0; col<=(NCols-1)/2; col++)
  {
    freq_values[col] = float(col)*DataResolution/(float(NCols));
    //cout << "Col: " << col << endl;
    lastcol = col;
  }
  for(int col = (NCols)/2; col > 0; col--)
  {
    lastcol++;
    freq_values[lastcol] = -float(col)*DataResolution/(float(NCols));
    //cout << "Col is: " << lastcol << "and selection is: " << -col << endl;
  }

  return freq_values;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns an array2d that contains entries scaled by 1/f^beta where
// beta is the fractal scaling
// This is used in generation of psuedo-fractal surfaces using the
// Fourier synthesis method of fractal generation
// (e.g., http://bringerp.free.fr/Files/Captain%20Blood/Saupe87d.pdf p.105
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D<float> LSDRasterSpectral::get_frequency_scaling_array(float beta)
{
  Array2D<float> freq_scaling_array(NRows,NCols);
  float radial_freq;

  // get the frequencies
  vector<float> row_freqs = get_row_direction_frequencies_unshifted();
  vector<float> col_freqs = get_col_direction_frequencies_unshifted();

  // get the frequency scaling
  for(int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      radial_freq = sqrt( row_freqs[row]*row_freqs[row]+col_freqs[col]*col_freqs[col]);

      if (radial_freq == 0)
      {
        freq_scaling_array[row][col] = 0;
      }
      else
      {
        freq_scaling_array[row][col] = 1/pow(radial_freq,beta);
      }
    }
  }

  return freq_scaling_array;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This creates a fractal surface using the spectral method
// The method works as follows:
//  1) Generate a random surface
//  2) Perform DFT on this random surface
//  3) Scale the tranform (both real and imaginary parts) by 1/f^beta
//  4) Perform the inverse DFT.
//
//  This results in a pseudo fractal surface that can be used in comarison
//  with real topography
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterSpectral::generate_fractal_surface_spectral_method(float beta)
{
  // first generate a random field
  float range = 1.0;
  rewrite_with_random_values(range);

  // now get the frequency scaling
  Array2D<float> freq_scaling = get_frequency_scaling_array(beta);

  // now get the real and imaginary DFT arrays
  // first set up the input and output arrays
  Array2D<float> InputArray = RasterData.copy();
  Array2D<float> OutputArrayReal = RasterData.copy();
  Array2D<float> OutputArrayImaginary = RasterData.copy();
  int transform_direction = -1;   // this means it will be a forward transform

  // perform the fourier annalysis
  dfftw2D_fwd(InputArray, OutputArrayReal, OutputArrayImaginary,transform_direction);
  cout << "Performed DFT! " << endl;

  // now scale the DFT. We replace values in the output arrays
  for(int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      OutputArrayReal[row][col] = OutputArrayReal[row][col]*freq_scaling[row][col];
      OutputArrayImaginary[row][col] = OutputArrayImaginary[row][col]*freq_scaling[row][col];
    }
  }

  // now perform the inverse transform
  // this overwrites the InputArray
  transform_direction = 1;
  dfftw2D_inv(OutputArrayReal, OutputArrayImaginary,
  	              InputArray, transform_direction);

  RasterData = InputArray.copy();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This creates a fractal surface using the spectral method
// The method works as follows:
//  1) Generate a random surface
//  2) Perform DFT on this random surface
//  3) Scale the tranform (both real and imaginary parts) by 1/f^beta
//  4) Perform the inverse DFT.
//
//  This results in a pseudo fractal surface that can be used in comarison
//  with real topography
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDRasterSpectral::generate_fractal_surface_spectral_method(float beta, float desired_relief)
{
  // first generate a random field
  float range = 1.0;
  rewrite_with_random_values(range);

  // now get the frequency scaling
  Array2D<float> freq_scaling = get_frequency_scaling_array(beta);

  // now get the real and imaginary DFT arrays
  // first set up the input and output arrays
  Array2D<float> InputArray = RasterData.copy();
  Array2D<float> OutputArrayReal = RasterData.copy();
  Array2D<float> OutputArrayImaginary = RasterData.copy();
  int transform_direction = -1;   // this means it will be a forward transform

  // perform the fourier annalysis
  dfftw2D_fwd(InputArray, OutputArrayReal, OutputArrayImaginary,transform_direction);
  cout << "Performed DFT! " << endl;

  // now scale the DFT. We replace values in the output arrays
  for(int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      OutputArrayReal[row][col] = OutputArrayReal[row][col]*freq_scaling[row][col];
      OutputArrayImaginary[row][col] = OutputArrayImaginary[row][col]*freq_scaling[row][col];
    }
  }

  // now perform the inverse transform
  // this overwrites the InputArray
  transform_direction = 1;
  dfftw2D_inv(OutputArrayReal, OutputArrayImaginary,
                InputArray, transform_direction);
                
                
                
  // now scale the data
  // now scale the DFT. We replace values in the output arrays
  float MaxElev = -999999999;
  float MinElev = 999999999;
  for(int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (InputArray[row][col] > MaxElev)
      {
        MaxElev = InputArray[row][col];
      }
      if (InputArray[row][col] < MinElev)
      {
        MinElev = InputArray[row][col];
      }
    }
  }
  
  float rel = MaxElev-MinElev;
  for(int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      InputArray[row][col] =  ((InputArray[row][col]-MinElev)/rel)*desired_relief;
    }
  }

  RasterData = InputArray.copy();
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//		FFFFFF  oOOo  UU  UU RRRR   IIIIII EEEEEE RRRR
//		FF     oO  Oo UU  UU RR RR    II   EE     RR RR
//		FFFF   OO  OO UU  UU RRRR     II   EEEE   RRRR
//		FF     oO  Oo UU  UU RR RR    II   EE     RR RR
//		FF      oOOo   uUUu  RR  RR IIIIII EEEEEE RR  RR
//
//     AAAA  NN    NN  AAAA  LL   YY    YY  sSSSs IIIIII  sSSSs
//    AA  AA NNNN  NN AA  AA LL    YY  YY  SS       II   SS
//    AAAAAA NN NN NN AAAAAA LL     YYYY    sSSs    II    sSSs
//    AA  AA NN  NNNN AA  AA LL      YY        SS   II       SS
//    AA  AA NN   NNN AA  AA LLLLLL  YY    SSSSs  IIIIII SSSSs
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// FAST FOURIER TRANSFORM MODULE
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Computes both the forward and inverse fast fourier transforms of a 2D
// discrete dataset.
// FOR FORWARD TRANSFORM:
//    - InputArray = zeta_padded (padded DEM)
//    - transform_direction = -1
//    - OutputArray = 2D spectrum
void LSDRasterSpectral::dfftw2D_fwd(Array2D<float>& InputArray, Array2D<float>& OutputArrayReal,
                                  Array2D<float>& OutputArrayImaginary, int transform_direction)
{
  //  cout << Ly << " " << Lx << endl;
  fftw_complex *input,*output;
  fftw_plan plan;

  // Declare one_dimensional contiguous arrays of dimension Ly*Lx
  input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Ly*Lx);
  output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Ly*Lx);

  // SET UP PLAN
  // -forwards, transform_direction==-1, -inverse, transform_direction==1
  if (transform_direction==-1)
  {
    //     cout << "  Running 2D discrete FORWARD fast fourier transform..." << endl;
    plan = fftw_plan_dft_2d(Ly,Lx,input,output,transform_direction,FFTW_MEASURE);
  }
  else
  {
    cout << "\nFATAL ERROR: for the tranform direction\n\t -1 = FORWARD \n\t" << endl;
		exit(EXIT_FAILURE);
  }

  // LOAD DATA INTO COMPLEX ARRAY FOR FFT IN ROW MAJOR ORDER
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      input[Lx*i+j][0] = InputArray[i][j];
      input[Lx*i+j][1] = 0;
    }
  }

  // EXECUTE PLAN
  fftw_execute(plan);

  // RETRIEVE OUTPUT - since data is real, we only need to extract real part of
  // the output.
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      OutputArrayReal[i][j] = output[Lx*i+j][0];
      OutputArrayImaginary[i][j] = output[Lx*i+j][1];
    }
  }

  // DEALLOCATE PLAN AND ARRAYS
  fftw_destroy_plan(plan);
  fftw_free(input);
  fftw_free(output);
}

//------------------------------------------------------------------------------
// FOR INVERSE TRANSFORM:
//    - InputArrays = Real and Imaginary components of 2D spectrum
//    - transform_direction = 1
//    - OutputArray = reconstructed DEM
void LSDRasterSpectral::dfftw2D_inv(Array2D<float>& InputArrayReal,
                          Array2D<float>& InputArrayImaginary,
                          Array2D<float>& OutputArray, int transform_direction)
{
  fftw_complex *input,*output;
  fftw_plan plan;

  // Declare one_dimensional contiguous arrays of dimension Ly*Lx
  input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Ly*Lx);
  output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Ly*Lx);

  // SET UP PLAN
  // -forwards => transform_direction==-1, -inverse => transform_direction==1
  if (transform_direction==1)
  {
//     cout << "  Running 2D discrete INVERSE fast fourier transform..." << endl;
    plan = fftw_plan_dft_2d(Ly,Lx,input,output,transform_direction,FFTW_MEASURE);
  }
  else
  {
    cout << "\nFATAL ERROR: for the tranform direction\n\t 1 = INVERSE \n\t" << endl;
		exit(EXIT_FAILURE);
  }

  // LOAD DATA INTO COMPLEX ARRAY FOR FFT IN ROW MAJOR ORDER
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      input[Lx*i+j][0] = InputArrayReal[i][j];
      input[Lx*i+j][1] = InputArrayImaginary[i][j];
    }
  }

  // EXECUTE PLAN
  fftw_execute(plan);

  // RETRIEVE OUTPUT ARRAY
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      OutputArray[i][j] = output[Lx*i+j][0];
    }
  }

  // DEALLOCATE PLAN AND ARRAYS
  fftw_destroy_plan(plan);
  fftw_free(input);
  fftw_free(output);
}

//------------------------------------------------------------------------------
// FOR INVERSE TRANSFORM: PASSING SINGLE COMPLEX ARRAY
//    - InputArray = Complex Array2D with real and imaginary components of 2D spectrum
//    - transform_direction = 1
//    - OutputArray = reconstructed DEM
//
// DAV - 22/10/2014  (Copied from Dave's version above and tweaked)
//------------------------------------------------------------------------------
void LSDRasterSpectral::dfftw2D_inv_complex(Array2D< complex< float > >& InputArrayComplex,
Array2D<float>& OutputArray, int transform_direction)
{
  /// Note using a float-version: "fftwf_complex" rather than "fftw_complex", which is type double.
  fftw_complex *input,*output;     
  fftw_plan plan;

  // Declare one_dimensional contiguous arrays of dimension Ly*Lx
  input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Ly*Lx);
  output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Ly*Lx);

  // SET UP PLAN
  // -forwards => transform_direction==-1, -inverse => transform_direction==1
  if (transform_direction==1)
  {
    cout << "  Running 2D discrete INVERSE fast fourier transform..." << endl;
    plan = fftw_plan_dft_2d(Ly,Lx,input,output,transform_direction,FFTW_MEASURE);
  }
  else
  {
    cout << "\nFATAL ERROR: for the tranform direction\n\t 1 = INVERSE \n\t" << endl;
    exit(EXIT_FAILURE);
  }

  // LOAD DATA INTO COMPLEX ARRAY FOR FFT IN ROW MAJOR ORDER
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      // In this version, we have passed a full complex array, so at this
      // point we then access the real and imaginary parts, so as not to change too much code
      input[Lx*i+j][0] = InputArrayComplex[i][j].real(); 	
      input[Lx*i+j][1] = InputArrayComplex[i][j].imag();
    }
  }

  // EXECUTE PLAN
  fftw_execute(plan);

  // RETRIEVE OUTPUT ARRAY
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      OutputArray[i][j] = output[Lx*i+j][0];
    }
  }

  // DEALLOCATE PLAN AND ARRAYS
  fftw_destroy_plan(plan);
  fftw_free(input);
  fftw_free(output);
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// DETREND DATA MODULE
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// FIT PLANE BY LEAST SQUARES REGRESSION AND USE COEFFICIENTS TO DETERMINE
// LOCAL SLOPE ax + by + c = z
// Have N simultaneous linear equations, and N unknowns.
// => b = Ax, where x is a 1xN array containing the coefficients we need for
// surface fitting.
// A is constructed using different combinations of x and y, thus we only need
// to compute this once, since the window size does not change.
// For 1st order surface fitting, there are 3 coefficients, therefore A is a
// 3x3 matrix
// Module kicks out detrended array, and an array with the trend plane
void LSDRasterSpectral::detrend2D(Array2D<float>& zeta, Array2D<float>& zeta_detrend, Array2D<float>& trend_plane)
{
//   cout << "  Detrending the DEM by fitting a planar surface..." << endl;
  Array2D<float> A(3,3,0.0);
  Array1D<float> bb(3,0.0);
  Array1D<float> coeffs(3);
  
  for (int i=0; i<NRows; ++i)
  {      
    for (int j=0; j<NCols; ++j)
    {
      if(zeta[i][j] != NoDataValue)
      {
        float x = j;
        float y = i;
        // Generate matrix A
        A[0][0] += pow(x,2);
        A[0][1] += x*y;
        A[0][2] += x;
        A[1][0] += y*x;
        A[1][1] += pow(y,2);
        A[1][2] += y;
        A[2][0] += x;
        A[2][1] += y;
        A[2][2] += 1;

        // Generate vector bb
        bb[0] += zeta[i][j]*x;
        bb[1] += zeta[i][j]*y;
        bb[2] += zeta[i][j];
      }
    }
  }

  // Solve matrix equations using LU decomposition using the TNT JAMA package:
  // A.coefs = b, where coefs is the coefficients vector.
  LU<float> sol_A(A);  // Create LU object
  coeffs = sol_A.solve(bb);
  float a_plane = coeffs[0];
  float b_plane = coeffs[1];
  float c_plane = coeffs[2];

  // Create detrended surface
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      float x = j;
      float y = i;
      trend_plane[i][j] = a_plane*x + b_plane*y + c_plane;

      if(zeta[i][j] != NoDataValue)
      {
        zeta_detrend[i][j] = zeta[i][j] - trend_plane[i][j];
      }
      else // Set NoDataValues as 0 on detrended surface
      {
        zeta_detrend[i][j] = 0;
      }
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// WINDOW MODULE
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Use 2D elliptical window on data matrix, to reduce spectral leakage and
// retain good frequency resolution. Return windowed data, the window itself and
// also the summed square of the weighting coefficients, WSS.
// Currently implemented are the unity window function (option 0), which is
// essentially no windowing and is the default, the Hann window(option = 1) and 
// Hamming window (window_option = 2). Another option would be to use a 2D Welch
// window, but this is yet to be coded up.
void LSDRasterSpectral::window_data(Array2D<float>& input, Array2D<float>& output, Array2D<float>& window, int window_option)
{
  if(window_option > 2)
  {
    cout << "window_option not available, default is a unity window" << endl;
    window_option =0;
  }
  // UNITY WINDOW
  if(window_option == 0)
  {
    output=input.copy();
    Array2D<float> window_array(NRows,NCols,1.0);
    window = window_array.copy(); 
    WSS=float(NRows*NCols); 
  }
  else
  {
    Array2D<float> output_array(NRows,NCols,0.0);
    Array2D<float> window_array(NRows,NCols,0.0);
    float PI = 3.14159265;
    float ny = float(NRows);
    float nx = float(NCols);
    // Get matrix coordinates of centroid of matrix
    float a = (nx-1)/2;
    float b = (ny-1)/2;
    WSS=0;
    
    // HANN WINDOW
    if(window_option == 1)
    {
      Array2D<float> r_prime_matrix(NRows,NCols,0.0);
      Array2D<float> theta_matrix(NRows,NCols,0.0);
      float r; // radial polar coordinate
      float theta; // angular polar coordinate
      float rprime;
      float HannCoefficient = 0;
      for(int i = 0; i < NRows; ++i)
      {
        for(int j = 0; j < NCols; ++j)
        {
          float x = j;
          float y = i;
          if(x == a)
          {
            theta = (PI/2);
          }
          else
          {
            theta = atan2((y - b),(x - a));
          }
          r = sqrt((y - b)*(y - b) + (x - a)*(x - a)); // distance from centre to this point
          
          // distance from centre to edge of ellipse for this particular theta
          rprime = sqrt((a*a)*(b*b)/(b*b*(cos(theta)*cos(theta)) + a*a*(sin(theta)*sin(theta)))); 
          if(r < rprime)
          {
            HannCoefficient = 0.5 * (1 + cos(PI * r/rprime));
            window_array[i][j] = HannCoefficient;
            WSS += HannCoefficient*HannCoefficient;
            output_array[i][j] = input[i][j] * HannCoefficient;
          }
        }
      }
    }
    // HAMMING WINDOW
    if(window_option==2)
    {
      float alpha = 0.54;
      float beta = 1-alpha;
      Array2D<float> r_prime_matrix(NRows,NCols,0.0);
      Array2D<float> theta_matrix(NRows,NCols,0.0);
      float r; // radial polar coordinate
      float theta; // angular polar coordinate
      float rprime;
      float HammingCoefficient = 0;
      for(int i = 0; i < NRows; ++i)
      {
        for(int j = 0; j < NCols; ++j)
        {
          float x = j;
          float y = i;
          if(x == a)
          {
            theta = (PI/2);
          }
          else
          {
            theta = atan2((y - b),(x - a));
          }
          r = sqrt((y - b)*(y - b) + (x - a)*(x - a)); // distance from centre to this point
          
          // distance from centre to edge of ellipse for this particular theta
          rprime = sqrt((a*a)*(b*b)/(b*b*(cos(theta)*cos(theta)) + a*a*(sin(theta)*sin(theta)))); 
          if(r < rprime)
          {
            HammingCoefficient = (alpha - beta*(PI * r/rprime));
            window_array[i][j] = HammingCoefficient;
            WSS += HammingCoefficient*HammingCoefficient;
            output_array[i][j] = input[i][j] * HammingCoefficient;
          }
        }
      }    
    }
//     // WELCH
//     if(window_option==3)
//     {
//       float r; // radial polar coordinate
//       float theta; // angular polar coordinate
//       float rprime;
//       float WelchCoefficient = 0;
//       for(int i = 0; i < NRows; ++i)
//       {
//         for(int j = 0; j < NCols; ++j)
//         {
//           float x = j;
//           float y = i;    
//           if(x == a)
//           {
//             theta = (PI/2);
//           }
//           else
//           {
//             theta = atan2((y - b),(x - a));
//           }
//           r = sqrt((y - b)*(y - b) + (x - a)*(x - a)); // distance from centre to this point
//           rprime = sqrt((a*a)*(b*b)/(b*b*(cos(theta)*cos(theta)) + a*a*(sin(theta)*sin(theta)))); // distance from centre to edge of ellipse for this particular theta
//           if(r < rprime)
//           {
//             WelchCoefficient = ;
//             window_array[i][j] = WelchCoefficient;
//             WSS += WelchCoefficient*WelchCoefficient;
//             output_array[i][j] = input[i][j] * WelchCoefficient;
//           }
//         }
//       }    
//     }
    output = output_array.copy();
    window = window_array.copy();
  }
  

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// HANN WINDOW MODULE
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Use 2D elliptical Hann (raised cosine) window on data matrix, to reduce
// spectral leakage and retain good frequency resolution.
// Return windowed data, the Hann window and also the summed square of the
// weighting coefficients, WSS.
// Another option would be to use a 2D Welch window, but functionality is very
// similar.
void LSDRasterSpectral::window_data_Hann2D(Array2D<float>& zeta_detrend, Array2D<float>& zeta_Hann2D, Array2D<float>& Hann2D)
{
  float PI = 3.14159265;
  //   cout << "  Windowing DEM using an elliptical 2D Hann window..." << endl;
  float ny = NRows;
  float nx = NCols;
  // Get matrix coordinates of centroid of matrix
  float a = (nx-1)/2;
  float b = (ny-1)/2;
  // Set up data window

  Array2D<float> r_prime_matrix(NRows,NCols,0.0);
  Array2D<float> id(NRows,NCols,0.0);
  Array2D<float> theta_matrix(NRows,NCols,0.0);
  float r; // radial polar coordinate
  float theta; // angular polar coordinate
  float rprime;
  float HannCoefficient = 0;
  WSS=0;
  for(int i = 0; i < NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      float x = j;
      float y = i;
      if(x == a)
      {
        theta = (PI/2);
      }
      else
      {
        theta = atan2((y - b),(x - a));
      }
      r = sqrt((y - b)*(y - b) + (x - a)*(x - a)); // distance from centre to this point
      rprime = sqrt((a*a)*(b*b)/(b*b*(cos(theta)*cos(theta)) + a*a*(sin(theta)*sin(theta)))); // distance from centre to edge of ellipse for this particular theta
      if(r < rprime)
      {
        HannCoefficient = 0.5 * (1 + cos(PI * r/rprime));
        Hann2D[i][j] = HannCoefficient;
        WSS += HannCoefficient*HannCoefficient;
        zeta_Hann2D[i][j] = zeta_detrend[i][j] * HannCoefficient;
        id[i][j]=1;
      }
    }
  }

}
// HAMMING WINDOW MODULE
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Use 2D elliptical Hann (raised cosine) window on data matrix, to reduce
// spectral leakage and retain good frequency resolution.
// Return windowed data, the Hann window and also the summed square of the
// weighting coefficients, WSS.
// Another option would be to use a 2D Welch window, but functionality is very
// similar.
void LSDRasterSpectral::window_data_Hamming2D(Array2D<float>& zeta_detrend, Array2D<float>& zeta_Hamming2D, Array2D<float>& Hamming2D)
{
  float PI = 3.14159265;
//   cout << "  Windowing DEM using an elliptical 2D Hann window..." << endl;
  float ny = NRows;
  float nx = NCols;
  // Get matrix coordinates of centroid of matrix
  float a = (nx-1)/2;
  float b = (ny-1)/2;
  // Set up data window
  float alpha = 0.54;
  float beta = 1-alpha;
  Array2D<float> r_prime_matrix(NRows,NCols,0.0);
  Array2D<float> id(NRows,NCols,0.0);
  Array2D<float> theta_matrix(NRows,NCols,0.0);
  float r; // radial polar coordinate
  float theta; // angular polar coordinate
  float rprime;
  float HammingCoefficient = 0;
  WSS=0;
  for(int i = 0; i < NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      float x = j;
      float y = i;
      if(x == a)
      {
        theta = (PI/2);
      }
      else
      {
        theta = atan2((y - b),(x - a));
      }
      r = sqrt((y - b)*(y - b) + (x - a)*(x - a)); // distance from centre to this point
      rprime = sqrt((a*a)*(b*b)/(b*b*(cos(theta)*cos(theta)) + a*a*(sin(theta)*sin(theta)))); // distance from centre to edge of ellipse for this particular theta
      if(r < rprime)
      {
        HammingCoefficient = (alpha - beta*(PI * r/rprime));
        Hamming2D[i][j] = HammingCoefficient;
        WSS += HammingCoefficient*HammingCoefficient;
        zeta_Hamming2D[i][j] = zeta_detrend[i][j] * HammingCoefficient;
        id[i][j]=1;
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// SHIFT ORIGIN OF SPECTRUM IN FOURIER DOMAIN
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// "The output of most algorithms that compute the DFT must be rearranged to
// place the zero wavenumber element near the center of the array. Provided Nx
// and Ny are even, dividing the output array into four equal quadrants and
// exchanging the nonadjacent quadrants will place the zero wavenumber element
// at the position (Nx/2, Ny/2) in the new array."  (Perron et al., 2008)

void LSDRasterSpectral::shift_spectrum(Array2D<float>& spectrum_real,  Array2D<float>& spectrum_imaginary, Array2D<float>& spectrum_real_shift, Array2D<float>& spectrum_imaginary_shift)
{
  vector<int> shifted_row_indices, shifted_col_indices;
  int p2y =(Ly+1)/2; 
  int p2x =(Lx+1)/2; 
  
  for(int i = p2y; i < Ly; ++i)  shifted_row_indices.push_back(i);
  for(int i = 0; i < p2y; ++i)   shifted_row_indices.push_back(i);
  for(int i = p2x; i < Lx; ++i)  shifted_col_indices.push_back(i);
  for(int i = 0; i < p2x; ++i)   shifted_col_indices.push_back(i);
  
  for(int row = 0; row < Ly; ++row)
  {
    for(int col = 0; col < Lx; ++col)
    {
       spectrum_real_shift[shifted_row_indices[row]][shifted_col_indices[col]] = spectrum_real[row][col];
       spectrum_imaginary_shift[shifted_row_indices[row]][shifted_col_indices[col]] = spectrum_imaginary[row][col];
    }
  }
}
void LSDRasterSpectral::shift_spectrum(Array2D<float>& spectrum_real,  Array2D<float>& spectrum_imaginary)
{
  Array2D<float> spectrum_real_shift(Ly,Lx,0.0);
  Array2D<float> spectrum_imaginary_shift(Ly,Lx,0.0);
  shift_spectrum(spectrum_real,  spectrum_imaginary, spectrum_real_shift, spectrum_imaginary_shift);
  spectrum_real=spectrum_real_shift.copy();
  spectrum_imaginary=spectrum_imaginary_shift.copy();
}

//------------------------------------------------------------------------------
// DE-SHIFT ORIGIN OF SPECTRUM
// Inverse process of above to return filtered spectrum to original format
// required for the inverse fourier transform algorithm.
void LSDRasterSpectral::shift_spectrum_inv(Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary, Array2D<float>& FilteredSpectrumReal_deshift, Array2D<float>& FilteredSpectrumImaginary_deshift)
{
  vector<int> i_shifted_row_indices, i_shifted_col_indices;
  int p2y =Ly-(Ly+1)/2; 
  int p2x =Lx-(Lx+1)/2;
  
  for(int i = p2y; i < Ly; ++i)  i_shifted_row_indices.push_back(i);
  for(int i = 0; i < p2y; ++i)   i_shifted_row_indices.push_back(i);
  for(int i = p2x; i < Lx; ++i)  i_shifted_col_indices.push_back(i);
  for(int i = 0; i < p2x; ++i)   i_shifted_col_indices.push_back(i);
  
  for(int row = 0; row < Ly; ++row)
  {
    for(int col = 0; col < Lx; ++col)
    {
       FilteredSpectrumReal_deshift[i_shifted_row_indices[row]][i_shifted_col_indices[col]] = FilteredSpectrumReal[row][col];
       FilteredSpectrumImaginary_deshift[i_shifted_row_indices[row]][i_shifted_col_indices[col]] = FilteredSpectrumImaginary[row][col];
    }
  }
}
void LSDRasterSpectral::shift_spectrum_inv(Array2D<float>& spectrum_real,  Array2D<float>& spectrum_imaginary)
{
  Array2D<float> spectrum_real_deshift(Ly,Lx,0.0);
  Array2D<float> spectrum_imaginary_deshift(Ly,Lx,0.0);
  shift_spectrum_inv(spectrum_real,  spectrum_imaginary, spectrum_real_deshift, spectrum_imaginary_deshift);
  spectrum_real=spectrum_real_deshift.copy();
  spectrum_imaginary=spectrum_imaginary_deshift.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// CALCULATE THE DFT PERIODOGRAM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Multiply fourier analysis output by complex conjugate and normalise.
// Note that for complex number z=x+iy, z*=x-iy, z.z* = x^2 + y^2
// Returns 2D PSD as only output
void LSDRasterSpectral::calculate_2D_PSD(Array2D<float>& spectrum_real_shift, Array2D<float>& spectrum_imaginary_shift)
{
  Array2D<float> temp(Ly,Lx,0.0);
  P_DFT = temp.copy();
  float fLx = float(Lx);
  float fLy = float(Ly);

  for (int i=0; i<Ly; ++i)
  {
    for (int j=0; j<Lx; ++j)
    {
      P_DFT[i][j] = fabs(pow(spectrum_real_shift[i][j],2) + pow(spectrum_imaginary_shift[i][j],2))/(fLy*fLx*WSS);
    }
  }
}
void LSDRasterSpectral::scale_spectrum(Array2D<float> SpectrumReal, Array2D<float> SpectrumIm,  float beta)
{
  float f,x,y;
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < Lx; ++j)
    {
      x = float(j);
      y = float(i);
      f = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx); // Radial Frequency
      if (f == 0)
      {
        SpectrumReal[i][j] = 0;
        SpectrumIm[i][j] = 0;
      }
      else
      {
        SpectrumReal[i][j] = SpectrumReal[i][j]/pow(f,beta);
        SpectrumIm[i][j] = SpectrumIm[i][j]/pow(f,beta);    
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// GET RADIAL POWER SPECTRUM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Collapse 2D PSD into a radial PSD
void LSDRasterSpectral::calculate_radial_PSD()
{
  // CALCULATE FREQUENCY INCREMENTS - for generation of power spectrum
  // Frequency goes from zero to 1/(2*resolution), the Nyquist frequency in
  // NRows_padded/2 increments.
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  vector<float> RadialFrequencyRaw(Ly*(Lx/2+1),0.0); // This is the distance from the origin in Frequency space.  Note that half of spectrum is redundant, since the fourier transform of a real dataset is symmetric, with a degeneracy of two.
  vector<float> RadialPSDRaw(Ly*(Lx/2+1),0.0);
//   float NyquistFreq = 1/(2*DataResolution);
  float RadialFreq;
  int count = 0;
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < (Lx/2+1); ++j)
    {
      float x = float(j);
      float y = float(i);
      float fLy = float(Ly);
      float fLx = float(Lx);
      RadialFreq = sqrt((y - (fLy/2))*(y - (fLy/2))*dfy*dfy + (x - (fLx/2))*(x - (fLx/2))*dfx*dfx); // distance from centre to this point. Converting position in frequency into an absolute frequency
      if (RadialFreq <= NyquistFreq)  // Ignore radial frequencies greater than the Nyquist frequency as these are aliased
        {
          RadialFrequencyRaw[count] = RadialFreq;
          RadialPSDRaw[count] = 2*P_DFT[i][j];   // Due to degeneracy
          ++count;
        }
    }
  }
  // Sort radial frequency
  vector<size_t> index_map;
  matlab_float_sort(RadialFrequencyRaw,RadialFrequencyRaw,index_map);
  // Reorder amplitudes to match sorted frequencies
  matlab_float_reorder(RadialPSDRaw,index_map,RadialPSDRaw);

//   // Get number of discrete radial frequencies
//   int n_freqs = 0;
//   for (int i=0; i<(Ly*(Lx/2+1)); ++i)
//   {
//     if (RadialFrequencyRaw[i] != RadialFrequencyRaw[i+1])
//     {
//       ++n_freqs;
//     }
//   }
// 
//   // Convert to spatially averaged spectrum
//   cout << "  Converting to radially averaged PSD..." << endl;
//   vector<float> RadialFrequency_grouped(n_freqs,0.0); // This is the distance from the origin in Frequency space
//   vector<float> RadialPSD_average(n_freqs,0.0);       // This will ultimately contain the radially averaged PSD
// 
//   int n_occurences = 0;           // This will keep track of the number of occurences of each radial frequency
//   int pointer = 0;
//   for (int i=0; i<(Ly*(Lx/2+1)); ++i)
//   {
//     RadialFrequency_grouped[pointer] = RadialFrequencyRaw[i];
//     RadialPSD_average[pointer] += RadialPSDRaw[i];
//     ++n_occurences;
// 
//     if (RadialFrequencyRaw[i] != RadialFrequencyRaw[i+1])
//     {
//       RadialPSD_average[pointer] = RadialPSD_average[pointer]/n_occurences;
//       // increment pointer and reset n_occurences
//       ++pointer;
//       n_occurences = 0;
//     }
//   }
//   RadiallyAveragedPSD = RadialPSD_average;
//   RadialFrequency = RadialFrequency_grouped;
  RadialPSD=RadialPSDRaw;
  RadialFrequency=RadialFrequencyRaw;   
//   // Copy across to output vectors
//   RadialFrequency_output = RadialFrequency_grouped;
//   RadialPSD_output = RadialPSD_average;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// MAIN FUNCTIONS USING SPECTRAL ANALYSIS
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// COMPUTE DISCRETE FAST FOURIER TRANSFORM OF A REAL, 2-DIMENSIONAL DATASET.
// Computes the 2D and radial power spectra of a 2D array.
// Input arguement is the width of the logarithmically spaced bins. For
// topography, suggest this is 0.1
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterSpectral::fftw2D_spectral_analysis(char* file_id, float LogBinWidth)
{
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DETREND DATA
  // FIT PLANE BY LEAST SQUARES REGRESSION AND USE COEFFICIENTS TO DETERMINE
  // LOCAL SLOPE ax + by + c = z
  Array2D<float> zeta_detrend(NRows,NCols);
  Array2D<float> trend_plane(NRows,NCols);
  detrend2D(RasterData, zeta_detrend, trend_plane);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // USE ELLIPTICAL 2D HANN (raised cosine) WINDOW ON ZETA MATRIX.
  // RETURN WINDOWED DATA AND ALSO THE SUMMED SQUARE OF THE WEIGHTING
  // COEFFICIENTS.
  Array2D<float> Hann2D(NRows,NCols,0.0);
  Array2D<float> zeta_Hann2D(NRows,NCols,0.0);
  WSS = 0; // summed square of weighting coefficients
  window_data_Hann2D(zeta_detrend, zeta_Hann2D, Hann2D);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // 2D DISCRETE FAST FOURIER TRANSFORM
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // PAD DATA WITH ZEROS TO A POWER OF TWO (facilitates FFT)

  Array2D<float> zeta_padded(Ly,Lx);
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      if (i<NRows && j<NCols)
      {
        zeta_padded[i][j] = zeta_Hann2D[i][j];
      }
      else
      {
        zeta_padded[i][j]=0;
      }
    }
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DO 2D FORWARD FAST FOURIER TRANSFORM
  int transform_direction = -1;
  Array2D<float> SpectrumReal(Ly,Lx);
  Array2D<float> SpectrumImaginary(Ly,Lx);
  dfftw2D_fwd(zeta_padded, SpectrumReal, SpectrumImaginary, transform_direction);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // REARRANGE SPECTRUM SO THAT ORIGIN IS AT THE CENTRE
  Array2D<float> SpectrumReal_shift(Ly,Lx);
  Array2D<float> SpectrumImaginary_shift(Ly,Lx);
  shift_spectrum(SpectrumReal, SpectrumImaginary, SpectrumReal_shift, SpectrumImaginary_shift);

  // CALCULATE THE DFT PERIODOGRAM
  // Multiply output by complex conjugate and normalise.
  // Note that for complex number z=x+iy, z*=x-iy, z.z* = x^2 + y^2
  calculate_2D_PSD(SpectrumReal_shift, SpectrumImaginary_shift);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // GET RADIAL POWER SPECTRUM
  // For forward transform, return the spectral power of the topography both
  // in a 2D array, and also as a one dimensional array of radial frequency
  calculate_radial_PSD();

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // BIN POWER SPECTRUM INTO LOGARITMICALLY SPACED BINS OF RADIAL FREQUENCY TO
  // GET MODEL "SIGNAL" FOR WIENER FILTER
  cout << "  Binning radial PSD into logarithmically spaced bins..." << endl;
  // Initiate output vectors
  vector<float> Bin_MeanRadialFreq;
  vector<float> Bin_RadialPSD;
  vector<float> BinMidpoints;
  vector<float> StandardDeviationRadialFreq;
  vector<float> StandardDeviationRadialPSD;
//   // Execute log binning
//   log_bin_data(RadialFrequency, RadiallyAveragedPSD, LogBinWidth, Bin_MeanRadialFreq, Bin_RadialPSD, BinMidpoints,
//                  StandardDeviationRadialFreq, StandardDeviationRadialPSD, int(NoDataValue)); 
  log_bin_data(RadialFrequency, RadialPSD, LogBinWidth, Bin_MeanRadialFreq, Bin_RadialPSD, BinMidpoints,
               StandardDeviationRadialFreq, StandardDeviationRadialPSD, int(NoDataValue));

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // WRITE OUTPUTS TO FILE
  int len = strlen(file_id);
  cout << "  Writing output files..." << endl;

  // 2D PSD
  char *PSD_file = new char[len+6];
  strcpy(PSD_file,file_id);
  strcat(PSD_file,"_P_DFT");
  
  LSDRaster PowerSpectrum(Ly,Lx,-(Lx/2),(Lx/2-1),DataResolution,NoDataValue,P_DFT);
  PowerSpectrum.write_raster(PSD_file,"flt");
  //----------------------------------------------------------------------------

	ofstream ofs;
  //----------------------------------------------------------------------------
  // Radially averaged PSD

  char *RadialPSD_file = new char[len+14];
  strcpy(RadialPSD_file,file_id);
  strcat(RadialPSD_file,"_RadialPSD.txt");

  ofs.open(RadialPSD_file);
  if( ofs.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << RadialPSD_file << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "Freq Wavelength PSD Model_PSD Model_noise\n";
  for(int i=0; i < int(RadialFrequency.size()); ++i)
  {
//     ofs << RadialFrequency[i] << " " << 1/RadialFrequency[i] << " " << RadiallyAveragedPSD[i] << " \n";
    ofs << RadialFrequency[i] << " " << 1/RadialFrequency[i] << " " << RadialPSD[i] << " \n";
  }
  ofs.close();
  //----------------------------------------------------------------------------
  // Binned averaged PSD
  char *RadialPSD_binned_file = new char[len+21];
  strcpy(RadialPSD_binned_file,file_id);
  strcat(RadialPSD_binned_file,"_RadialPSD_binned.txt");

  ofs.open(RadialPSD_binned_file);
  if( ofs.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << RadialPSD_binned_file << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "Freq Wavelength PSD Sigma Model Noise\n";
  for(int i=0; i < int(Bin_MeanRadialFreq.size()); ++i)
  {
    ofs << Bin_MeanRadialFreq[i] << " " << 1/Bin_MeanRadialFreq[i] << " " << Bin_RadialPSD[i] << " " << StandardDeviationRadialPSD[i] << " \n";
  }
  ofs.close();

  //----------------------------------------------------------------------------
  cout << "  DONE!" << endl;
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//		 sSSSs PPPPP  EEEEEE   CCCC TTTTTT RRRR    AAAA  LL
//		SS     PP  PP EE      CC      TT   RR RR  AA  AA LL
//		 sSSs  PPPPP  EEEE   CC       TT   RRRR   AAAAAA LL
//		    SS PP     EE      CC      TT   RR RR  AA  AA LL
//		SSSSs  PP     EEEEEE   CCCC   TT   RR  RR AA  AA LLLLLL
//
//    FFFFFF IIIIII LL   TTTTTT EEEEEE RRRR    sSSSs
//    FF       II   LL     TT   EE     RR RR  SS
//    FFFF     II   LL     TT   EEEE   RRRR    sSSs
//    FF       II   LL     TT   EE     RR RR      SS
//    FF     IIIIII LLLLLL TT   EEEEEE RR  RR SSSSs
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// FILTER WEIGHTS
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// BANDPASS FILTER
// Filter array to band between frequency bands f1 and f2.  The bandpass filter
// is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
void LSDRasterSpectral::bandpass_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
                                        Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
                                        float f1, float f2)
{
  // CALCULATE FREQUENCY INCREMENTS - for generation of power spectrum
  // Frequency goes from zero to 1/(2*resolution), the Nyquist frequency in
  // NRows_padded/2 increments.
//   float dfx = 1/(DataResolution*Lx);
//   float dfy = 1/(DataResolution*Ly);
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cout << "  Gaussian bandpass filter between f1 = " << f1 << " and f2 = " << f2 << endl;
  float f; // radial frequency
  float weight; // Filter weight
  float sigma = sqrt((f2-f1)*(f2-f1))/6; // Standard Deviation of Gaussian filter
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < Lx; ++j)
    {
      float x = j;
      float y = i;
      // Converting position in frequency space into an absolute frequency
      f = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx);
      weight = exp(-(f - 0.5*(f1 + f2))*(f - 0.5*(f1 + f2))/(2*sigma*sigma));
      FilteredSpectrumReal[i][j] = weight*RawSpectrumReal[i][j];
      FilteredSpectrumImaginary[i][j] = weight*RawSpectrumImaginary[i][j];
    }
  }
}
//------------------------------------------------------------------------------
// LOWPASS FILTER
// Filter array to retain frequencies below f1.  The filter edge is a radial
// gaussian function with a SD of |f2-f1|/3.
void LSDRasterSpectral::lowpass_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
                                       Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
                                       float f1, float f2)
{
  // CALCULATE FREQUENCY INCREMENTS - for generation of power spectrum
  // Frequency goes from zero to 1/(2*resolution), the Nyquist frequency in
  // NRows_padded/2 increments.
  //----------------------------------------------------------------------------
  cout << "  Lowpass filter with edges controlled by radial Gaussian function between f1 = " << f1 << " and f2 = " << f2 << endl;
  float f; // radial frequency
  float weight; // Filter weight
  float sigma;   // Standard Deviation of Gaussian edge
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < Lx; ++j)
    {
      float x = j;
      float y = i;
      f = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx);
                     // distance from centre to this point. Converting position in frequency into an absolute frequency
      if (f < f1)
      {
        weight = 1;
      }
      else
      {
        if (f2 > f1)
        {
          sigma = sqrt((f2-f1)*(f2-f1))/3;
          weight = exp(-(f - f1)*(f-f1)/(2*sigma*sigma));
        }
        else
        {
          weight = 0; // this is for the case that f1 = f2 and essentially the weighting function acts with a hard edge at f=f1=f2.
        }
      }
      FilteredSpectrumReal[i][j] = weight*RawSpectrumReal[i][j];
      FilteredSpectrumImaginary[i][j] = weight*RawSpectrumImaginary[i][j];
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LOWPASS FILTER REMAINDER
// Filter array to retain frequencies below f1.  The filter edge is a radial
// gaussian function with a SD of |f2-f1|/3.
void LSDRasterSpectral::lowpass_filter_remainder(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
                                       Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
                                       float f1, float f2)
{
  // CALCULATE FREQUENCY INCREMENTS - for generation of power spectrum
  // Frequency goes from zero to 1/(2*resolution), the Nyquist frequency in
  // NRows_padded/2 increments.
  //----------------------------------------------------------------------------
  cout << "  Lowpass filter with edges controlled by radial Gaussian function between f1 = " << f1 << " and f2 = " << f2 << endl;
  float f; // radial frequency
  float weight; // Filter weight
  float sigma;   // Standard Deviation of Gaussian edge
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < Lx; ++j)
    {
      float x = j;
      float y = i;
      f = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx);
      		               // distance from centre to this point. Converting position in frequency into an absolute frequency
      if (f < f1)
      {
        weight = 0;
      }
      else
      {
        if (f2 > f1)
        {
          sigma = sqrt((f2-f1)*(f2-f1))/3;
          weight = 1-exp(-(f - f1)*(f-f1)/(2*sigma*sigma));
        }
        else
        {
          weight = 1; // this is for the case that f1 = f2 and essentially the weighting function acts with a hard edge at f=f1=f2.
        }
      }
      FilteredSpectrumReal[i][j] = weight*RawSpectrumReal[i][j];
      FilteredSpectrumImaginary[i][j] = weight*RawSpectrumImaginary[i][j];
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//------------------------------------------------------------------------------
// HIGHPASS FILTER
// Filter array to retain frequencies above f1.  The filter edge is a radial
// gaussian function with a SD of |f2-f1|/3.
void LSDRasterSpectral::highpass_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
                                        Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
                                        float f1, float f2)
{
  // CALCULATE FREQUENCY INCREMENTS - for generation of power spectrum
  // Frequency goes from zero to 1/(2*resolution), the Nyquist frequency in
  // NRows_padded/2 increments.
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cout << "    Highpass filter with edges controlled by radial Gaussian function between f1 = " << f1 << " and f2 = " << f2 << endl;
  float f; // radial frequency
  float weight; // Filter weight
  float sigma = sqrt((f2-f1)*(f2-f1))/3; // Standard Deviation of Gaussian filter
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < Lx; ++j)
    {
      float x = j;
      float y = i;
      f = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx); // distance from centre to this point. Converting position in frequency into an absolute frequency
      if (f > f2)
      {
        weight = 1;
      }
      else
      {
        if (f2 > f1)
        {
          sigma = sqrt((f2-f1)*(f2-f1))/3;
          weight = exp(-(f - f2)*(f-f2)/(2*sigma*sigma));
        }
        else
        {
          weight = 0; // this is for the case that f1 = f2 and essentially the weighting function acts with a hard edge at f=f1=f2.
        }
      }
      FilteredSpectrumReal[i][j] = weight*RawSpectrumReal[i][j];
      FilteredSpectrumImaginary[i][j] = weight*RawSpectrumImaginary[i][j];
    }
  }
}
//------------------------------------------------------------------------------
// HIGHPASS FILTER REMAINDER
// Filter weights 1-FHigh -> this gives the component of the spectrum filtered
// out by the highpass filter
void LSDRasterSpectral::highpass_filter_remainder(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
                                        Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
                                        float f1, float f2)
{
  // CALCULATE FREQUENCY INCREMENTS - for generation of power spectrum
  // Frequency goes from zero to 1/(2*resolution), the Nyquist frequency in
  // NRows_padded/2 increments.
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cout << "    Lowpass filter version 2, this gives residual topography to that produced from a highpass filter with edges controlled by radial Gaussian function between f1 = " << f1 << " and f2 = " << f2 << endl;
  float f; // radial frequency
  float weight; // Filter weight
  float sigma = sqrt((f2-f1)*(f2-f1))/3; // Standard Deviation of Gaussian filter
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < Lx; ++j)
    {
      float x = j;
      float y = i;
      f = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx); // distance from centre to this point. Converting position in frequency into an absolute frequency
      if (f > f2)
      {
        weight = 0;
      }
      else
      {
        if (f2 > f1)
        {
          sigma = sqrt((f2-f1)*(f2-f1))/3;
          weight = 1 - exp(-(f - f2)*(f-f2)/(2*sigma*sigma));
        }
        else
        {
          weight = 1; // this is for the case that f1 = f2 and essentially the weighting function acts with a hard edge at f=f1=f2.
        }
      }
      FilteredSpectrumReal[i][j] = weight*RawSpectrumReal[i][j];
      FilteredSpectrumImaginary[i][j] = weight*RawSpectrumImaginary[i][j];
    }
  }
}

//------------------------------------------------------------------------------
// WIENER FILTER
// The Wiener filter is a spectral filter that removes noise from an image or
// DEM.  Essentially, it works on the principle that the observed spectrum
// contains the superposition of the real signal and an additional noise signal,
// which we want to remove.  If we know, or can make a reasonable guess at the
// noise, N(f), and signal, S(f), parts of the spectrum then we can remove the
// noise using the filter:
//
//        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
//
// For topography; at long wavelengths the topographic signal obeys an
// approximate power law relationship between amplitude and frequency,
// decreasing as the frequency increases (and wavelength decreases).  Noise
// typically dominates the high frequency part of the spectrum.  Thus at high
// frequencies the spectrum is dominated by noise, and the filter weight goes to
// zero.  In contrast, at low frequencies, the signal dominates and the filter
// weight goes to 1.
//
void LSDRasterSpectral::wiener_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary, Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary)
{
  // CALCULATE FREQUENCY INCREMENTS - for generation of power spectrum
  // Frequency goes from zero to 1/(2*resolution), the Nyquist frequency in
  // NRows_padded/2 increments.
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // GET 2D POWER SPECTRUM
  calculate_2D_PSD(RawSpectrumReal, RawSpectrumImaginary);
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // GET RADIAL POWER SPECTRUM
  // For forward transform, return the spectral power of the topography both
  // in a 2D array, and also as a one dimensional array of radial frequency
  calculate_radial_PSD();
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // FIT POWER LAW TO SPECTRUM BETWEEN RANGE OF WAVELENGTHS 1000m - 100m (THE
  // RANGE EXPECTED TO FALL WITHIN WAVELENGTHS CONTROLLED BY RIDGE-VALLEY
  // TOPOGRAPHY)
  cout << "  Fitting power law to data..." << endl;
  vector<float> LogRadialFrequency;
  vector<float> LogRadialPSD;
  int n_freqs = RadialFrequency.size();
  float f_low = 0.001; // frequency at wavelength of 1000m
  float f_high = 0.01; // frequency at wavelength of 100m
  float m_model,logc_model,c_model;      // Coefficients of power law fit => logPSD = logc + m*log(freq) => PSD = c*freq^m
  for (int i = 0; i < n_freqs; ++i)
  {
    //cout << RadialFrequency[i] << endl;
    if(RadialFrequency[i] <= f_high && RadialFrequency[i] >= f_low)
    {
      LogRadialFrequency.push_back(log10(RadialFrequency[i]));
//       LogRadialPSD.push_back(log10(RadiallyAveragedPSD[i]));
      LogRadialPSD.push_back(log10(RadialPSD[i]));
    }
  }
  // Least squares regression
  vector<float> residuals;
  cout << LogRadialFrequency.size() << " " << LogRadialPSD.size() << " " <<  endl;
  vector<float> regression_results = simple_linear_regression(LogRadialFrequency, LogRadialPSD, residuals);
  cout << "done" << endl;
  m_model = regression_results[0];
  logc_model = regression_results[1];

  //linear_fit(LogRadialFrequency, LogRadialPSD, m_model, logc_model);
  c_model = pow(10,logc_model);

  // Extend relationship across entire frequency range to produce model spectrum
  vector<float> RadialPSD_model(n_freqs,0.0);
  for (int i = 0; i < n_freqs; ++i)
  {
    if(RadialFrequency[i] == 0)
    {
      RadialPSD_model[i] = 0;
    }
    else
    {
      RadialPSD_model[i] = c_model*pow(RadialFrequency[i],m_model);
    }
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // GET MEAN AMPLITUDE OF SPECTRUM CLOSE TO NYQUIST FREQUENCY AS ESTIMATE OF
  // WHITE NOISE SPECTRUM

  // MAKE INTO IF LOOP WITH CONDITION FOR WHITE NOISE VS MODEL NOISE

  // Note that this assumes that high frequency "noise" is random.  We might
  // expect that the high frequency content might actually be more structured
  // than that, particularly if the high frequency signal is being produced by
  // rock exposures, unfiltered vegetation, or even pesky gophers.
  cout << "  Estimating noise distribution as white noise... " << endl;

  // i) Basically going to use a simple high-pass filter, and then average the
  // amplitude to give noise.
  //vector<float> WhiteNoise(n_freqs,0.0);
  int n_freqs_noise = 0;
  float f_highpass = pow(10,-0.7);
  
  // Make sure that the highpass filter is not greater than the DEM resolution. 
  // If it is, set the highpass filter to 3x the DEM resolution
  float L_highpass = 1/f_highpass;
  if(L_highpass < 2*DataResolution)
  {
    L_highpass = 2*DataResolution+0.00001;
    f_highpass = 1/L_highpass;
  }
  
  float WhiteNoiseAmplitude = 0;
//   for (int i = 0; i < int(RadiallyAveragedPSD.size()); ++i)
  for (int i = 0; i < int(RadialPSD.size()); ++i)
  {
    if(RadialFrequency[i] >= f_highpass)
    {
//       WhiteNoiseAmplitude += RadiallyAveragedPSD[i];
      WhiteNoiseAmplitude += RadialPSD[i];
      ++n_freqs_noise;
    }
  }
  WhiteNoiseAmplitude = WhiteNoiseAmplitude/n_freqs_noise;
  //----------------------------------------------------------------------------
  // ii) Alternatively can model noise using a linear fit through the high
  // frequency part of the spectrum - in particular it would be interesting to
  // see if it approximates pink noise, which is ubiquitous to natural systems
  // with self-organised criticality.
  //float m_noise,logc_noise,c_noise;
  //vector<float> LogRadialFrequency_highpass;
  //vector<float> LogRadialPSD_highpass;
  //for (int i = 0; i < n_freqs; ++i)
  //{
  //  if(RadialFrequency[i] >= f_highpass)
  //  {
  //    LogRadialFrequency_highpass.push_back(log10(RadialFrequency[i]));
  //    LogRadialPSD_highpass.push_back(log10(RadialPSD[i]));
  //  }
  //}
  //// Least squares regression
  //regression_results = simple_linear_regression(LogRadialFrequency_highpass, LogRadialPSD_highpass,  residuals);
  //m_noise = regression_results[0];
  //logc_noise = regression_results[1];


  //linear_fit(LogRadialFrequency_highpass, LogRadialPSD_highpass, m_noise, logc_noise);


  //c_noise = pow(10,logc_noise);
  // Extend relationship across entire frequency range to produce model spectrum
  vector<float> Noise_model(n_freqs,0.0);
  for (int i = 0; i < n_freqs; ++i)
  {
    //Noise_model[i] = c_noise*pow(RadialFrequency[i],m_noise);
    Noise_model[i] = WhiteNoiseAmplitude;
  }
  //cout << "Modeled noise exponent = " << m_noise << endl;
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // WIENER FILTER
  // Determine Wiener Coefficients and apply to spectrum
  // WienerCoefficient = Signal/(Signal + Noise).  Basically acts as a lowpass
  // filter to remove noise from image.
  float model;
  float noise;
  float f; // radial frequency
  float WienerCoefficient; // Filter weight
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < Lx; ++j)
    {
      float x = j;
      float y = i;
      f = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx); // Radial Frequency
      model = c_model*pow(f,m_model);
      //noise = c_noise*pow(f,m_noise);
      noise = WhiteNoiseAmplitude;
      if (f == 0) WienerCoefficient = 1;
      else WienerCoefficient = model/(model+noise);
      FilteredSpectrumReal[i][j] = WienerCoefficient*RawSpectrumReal[i][j];
      FilteredSpectrumImaginary[i][j] = WienerCoefficient*RawSpectrumImaginary[i][j];
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// MAIN FUNCTIONS USING SPECTRAL FILTERS
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// FAST FOURIER TRANSFORM FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
//
// Note that FLow <= FHigh
//
// There are three types of filters depending on the intentions of the user
//
// BANDPASS FILTER (FilterType = 1)
// Filter array to band between frequency bands f1 and f2.  The bandpass filter
// is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
//
// LOWPASS FILTER (FilterType = 2)
// Filter array to retain frequencies below f1.  The filter edge is a radial
// gaussian function with a SD of |f2-f1|/3.  f1 is the frequency below which
// the filter starts to taper; f2 is the frequency at which the filter tapers to
// zero. If f1 = f2, the edge is effectively a step function.
// HIGHPASS FILTER (FilterType = 3)
//
// Filter array to retain frequencies above f2.  The filter edge is a radial
// gaussian function with a SD of |f2-f1|/3.  f2 is the frequency below which
// the filter starts to taper; f1 is the frequency at which the filter tapers to
// zero. If f1 = f2, the edge is effectively a step function.
//
// A second type of bandpass filter is possible by combining the highpass and
// lowpass filters.
//------------------------------------------------------------------------------
// David Milodowski, 10/12/2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDRasterSpectral::fftw2D_filter(int FilterType, float FLow, float FHigh)
{
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // for FORWARD TRANSFORM
  cout << "\n***fftw_2Dfilt_v1.1: spectral filtering of array***" << endl;
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DETREND DATA => DO NOT WINDOW!
  // FIT PLANE BY LEAST SQUARES REGRESSION AND USE COEFFICIENTS TO DETERMINE
  // LOCAL SLOPE ax + by + c = z
  Array2D<float> zeta_detrend(NRows,NCols);
  Array2D<float> trend_plane(NRows,NCols);
  detrend2D(RasterData, zeta_detrend, trend_plane);
  //float WSS = 1; // dataset is not windowed

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // 2D DISCRETE FAST FOURIER TRANSFORM
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // PAD DATA WITH ZEROS TO A POWER OF TWO (facilitates FFT)

  Array2D<float> zeta_padded(Ly,Lx);
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      if (i<NRows && j<NCols)
      {
        zeta_padded[i][j] = zeta_detrend[i][j];
      }
      else
      {
        zeta_padded[i][j]=0;
      }
    }
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DO 2D FORWARD FAST FOURIER TRANSFORM
  int transform_direction = -1;
  Array2D<float> spectrum_real(Ly,Lx);
  Array2D<float> spectrum_imaginary(Ly,Lx);
  dfftw2D_fwd(zeta_padded, spectrum_real, spectrum_imaginary, transform_direction);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // REARRANGE SPECTRUM SO THAT ORIGIN IS AT THE CENTRE
  Array2D<float> spectrum_real_shift(Ly,Lx);
  Array2D<float> spectrum_imaginary_shift(Ly,Lx);
  shift_spectrum(spectrum_real, spectrum_imaginary, spectrum_real_shift, spectrum_imaginary_shift);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // INVERSE TRANSFORM
  // For inverse tranform, take 2D power spectrum.  Filter for desired frequncy
  // band and return spectrally filtered topography
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // APPLY FILTER
  // First, remove the frequency ranges that are not wanted
  Array2D<float> FilteredSpectrumReal(Ly,Lx,0.0);
  Array2D<float> FilteredSpectrumImaginary(Ly,Lx,0.0);

  // SET FILTER PARAMETERS
  if (FilterType == 1)
  {
    bandpass_filter(spectrum_real_shift, spectrum_imaginary_shift, FilteredSpectrumReal, FilteredSpectrumImaginary,FLow, FHigh);
  }
  else if (FilterType == 2)
  {
    lowpass_filter(spectrum_real_shift, spectrum_imaginary_shift, FilteredSpectrumReal, FilteredSpectrumImaginary, FLow, FHigh);
  }
  else if (FilterType == 3)
  {
    highpass_filter(spectrum_real_shift, spectrum_imaginary_shift, FilteredSpectrumReal, FilteredSpectrumImaginary, FLow, FHigh);
  }  
  else if (FilterType == 4)
  {
    lowpass_filter_remainder(spectrum_real_shift, spectrum_imaginary_shift, FilteredSpectrumReal, FilteredSpectrumImaginary, FLow, FHigh);
  }
  else if (FilterType == 5)
  {
    highpass_filter_remainder(spectrum_real_shift, spectrum_imaginary_shift, FilteredSpectrumReal, FilteredSpectrumImaginary, FLow, FHigh);
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DE-SHIFT ORIGIN OF SPECTRUM
  // Return filtered spectrum to original format (de-shifted)
  Array2D<float> FilteredSpectrumReal_deshift(Ly,Lx);
  Array2D<float> FilteredSpectrumImaginary_deshift(Ly,Lx);
  shift_spectrum_inv(FilteredSpectrumReal, FilteredSpectrumImaginary, FilteredSpectrumReal_deshift, FilteredSpectrumImaginary_deshift);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DO 2D INVERSE FAST FOURIER TRANSFORM
  transform_direction = 1;
  Array2D<float> FilteredTopographyPadded(Ly,Lx);
  dfftw2D_inv(FilteredSpectrumReal_deshift, FilteredSpectrumImaginary_deshift, FilteredTopographyPadded, transform_direction);
  // Need to scale output by the number of pixels, and by the Hann window to
  // recover the topography, before adding the planar trend back to the dataset
  cout << "  Scaling output filtered topography..." << endl;
  Array2D<float> FilteredTopography(NRows,NCols,NoDataValue);
  for (int i=0; i < NRows; ++i)
  {
    for (int j=0; j < NCols; ++j)
    {
      if(RasterData[i][j]!=NoDataValue) FilteredTopography[i][j] = FilteredTopographyPadded[i][j]/(Lx*Ly) + trend_plane[i][j];
    }
  }
  // Remove areas initially classified as NoData in the original DEM (only really necessary for irregular shaped datasets)
  for (int i=0; i < NRows; ++i)
  {
    for (int j=0; j < NCols; ++j)
    {
      if(RasterData[i][j]==NoDataValue) FilteredTopography[i][j] = NoDataValue;
    }
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  LSDRaster FilteredTopographyRaster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FilteredTopography,GeoReferencingStrings);
  return FilteredTopographyRaster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// WIENER FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// The Wiener filter is a spectral filter that removes noise from an image or
// DEM.  Essentially, it works on the principle that the observed spectrum
// contains the superposition of the real signal and an additional noise signal,
// which we want to remove.  If we know, or can make a reasonable guess at the
// noise, N(f), and signal, S(f), parts of the spectrum then we can remove the
// noise using the filter:
//
//        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
//
// For topography; at long wavelengths the topographic signal obeys an
// approximate power law relationship between amplitude and frequency,
// decreasing as the frequency increases (and wavelength decreases).  Noise
// typically dominates the high frequency part of the spectrum.  Thus at high
// frequencies the spectrum is dominated by noise, and the filter weight goes to
// zero.  In contrast, at low frequencies, the signal dominates and the filter
// weight goes to 1.
//
// The optimal wiener filter is described in more detail in Numerical Recipes,
// 13.3, p149.
//
// The exact structure of the noise is worth thinking about.  White noise, which
// is random, has equal power across all wavelengths.  In the instance of
// topography, noise can be created by a whole range of sources, from rock
// exposure, to pit and mound topography, to unfiltered vegetation etc.  It is
// likely that these sources will not produce purely white noise, but rather
// will show an element of structure.  This program makes two assumptions about
// the noise: i) it dominates the signal at high frequencies (close to the
// Nquist frequency) and ii) we can reasonably model this using a linear fit in
// log-log space - i.e. it obeys some form of power law function between
// frequency and amplitude.  Note that if the noise in the signal is really
// white noise, then the power law function for the noise would simply have an
// exponent of zero.  I prefer this formulation because it permits the
// characterisation of the noise model without assuming that the noise has a
// particular structure (white noise, pink noise etc.)
//
//------------------------------------------------------------------------------
// David Milodowski, 10/12/2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDRasterSpectral::fftw2D_wiener()
{
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DETREND DATA => DO NOT WINDOW!
  // FIT PLANE BY LEAST SQUARES REGRESSION AND USE COEFFICIENTS TO DETERMINE
  // LOCAL SLOPE ax + by + c = z
  Array2D<float> zeta_detrend(NRows,NCols);
  Array2D<float> trend_plane(NRows,NCols);
  detrend2D(RasterData, zeta_detrend, trend_plane);
  WSS = 1; // dataset is not windowed

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // 2D DISCRETE FAST FOURIER TRANSFORM
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // PAD DATA WITH ZEROS TO A POWER OF TWO (facilitates FFT)
  Ly = int(pow(2,ceil(log(NRows)/log(2))));
  Lx = int(pow(2,ceil(log(NCols)/log(2))));

  Array2D<float> zeta_padded(Ly,Lx);
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      if (i<NRows && j<NCols)
      {
        zeta_padded[i][j] = zeta_detrend[i][j];
      }
      else
      {
        zeta_padded[i][j]=0;
      }
    }
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DO 2D FORWARD FAST FOURIER TRANSFORM
  int transform_direction = -1;
  Array2D<float> spectrum_real(Ly,Lx);
  Array2D<float> spectrum_imaginary(Ly,Lx);
  dfftw2D_fwd(zeta_padded, spectrum_real, spectrum_imaginary, transform_direction);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // REARRANGE SPECTRUM SO THAT ORIGIN IS AT THE CENTRE
  Array2D<float> spectrum_real_shift(Ly,Lx);
  Array2D<float> spectrum_imaginary_shift(Ly,Lx);
  shift_spectrum(spectrum_real, spectrum_imaginary, spectrum_real_shift, spectrum_imaginary_shift);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // INVERSE TRANSFORM
  // For inverse tranform, take 2D power spectrum.  Filter for desired frequncy
  // band and return spectrally filtered topography
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // APPLY FILTER
  // First, remove the frequency ranges that are not wanted
  Array2D<float> FilteredSpectrumReal(Ly,Lx,0.0);
  Array2D<float> FilteredSpectrumImaginary(Ly,Lx,0.0);
  // SET FILTER PARAMETERS
  wiener_filter(spectrum_real_shift, spectrum_imaginary_shift, FilteredSpectrumReal, FilteredSpectrumImaginary);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DE-SHIFT ORIGIN OF SPECTRUM
  // Return filtered spectrum to original format (de-shifted)
  Array2D<float> FilteredSpectrumReal_deshift(Ly,Lx);
  Array2D<float> FilteredSpectrumImaginary_deshift(Ly,Lx);
  shift_spectrum_inv(FilteredSpectrumReal, FilteredSpectrumImaginary, FilteredSpectrumReal_deshift, FilteredSpectrumImaginary_deshift);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DO 2D INVERSE FAST FOURIER TRANSFORM
  transform_direction = 1;
  Array2D<float> FilteredTopographyPadded(Ly,Lx);
  dfftw2D_inv(FilteredSpectrumReal_deshift, FilteredSpectrumImaginary_deshift, FilteredTopographyPadded, transform_direction);
  // Need to scale output by the number of pixels, and by the Hann window to
  // recover the topography, before adding the planar trend back to the dataset
  cout << "  Scaling output filtered topography..." << endl;
  Array2D<float> FilteredTopography(NRows,NCols,NoDataValue);
  for (int i=0; i < NRows; ++i)
  {
    for (int j=0; j < NCols; ++j)
    {
      if(RasterData[i][j]!=NoDataValue) FilteredTopography[i][j] = FilteredTopographyPadded[i][j]/(Lx*Ly) + trend_plane[i][j];
    }
  }
  // Remove areas initially classified as NoData in the original DEM (only really necessary for irregular shaped datasets)
  for (int i=0; i < NRows; ++i)
  {
    for (int j=0; j < NCols; ++j)
    {
      if(RasterData[i][j]==NoDataValue) FilteredTopography[i][j] = NoDataValue;
    }
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  LSDRaster FilteredTopographyRaster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FilteredTopography,GeoReferencingStrings);
  return FilteredTopographyRaster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// SPECTRUM NORMALISATION - following Perron et al., 2008
// Functions to estimate roll-over frequency for landscapes, and calculate the
// background spectrum, which can then be used to normalise landscape spectra to
// facilitate the identification of emergent lengthscales.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// find_rollover_frequency
void LSDRasterSpectral::find_rollover_frequency(float& rollover_frequency, float& rollover_beta,float& sub_rollover_beta, float log_bin_width)
{
  // First bin the data into logarithmically spaced bins (default is 50).
  vector<float> bin_midpoints,bin_mean_PSD,bin_mean_freq;
  vector< vector<float> > logbinned_PSD;
  cout << "\t find rollover frequency\t 1) bin data" << endl;
  log_bin_data(RadialPSD, RadialFrequency, log_bin_width, bin_midpoints, bin_mean_PSD, bin_mean_freq, logbinned_PSD);
  
  // Loop through NBins - 1 times, running linear regression on the binned 
  // frequencies in log-log space.  in each iteration, removing bins of lower
  // frequency.  Calculate the R-squared value and beta value
  float beta_max = 0;
  float r2_max = 0;
  float rollover_freq = 0;
  //float rollover_freq_r2 = 0;
  int NBins = bin_mean_PSD.size();
  vector<float> log_bin_mean_PSD,log_bin_mean_freq;
  vector<int> bin_index;
  cout << "\t\t\t\t\t 2) repeated linear regressions" << endl;
  
  // get log of each bin mean
  for(int i = 0; i < NBins; ++i)
  {
    if((bin_mean_PSD[i]!=NoDataValue)&&(bin_mean_freq[i]!=NoDataValue))
    {
      if(bin_mean_PSD[i]==0 || bin_mean_freq[i] == 0)
      {
        cout << "Ahhh I am about to take the log of zero! Segmentation fault imminent!!!" << endl;
      }
      log_bin_mean_PSD.push_back(log10(bin_mean_PSD[i]));
      log_bin_mean_freq.push_back(log10(bin_mean_freq[i]));
      bin_index.push_back(i);
    }
  }
  
  int cutoff = log_bin_mean_freq.size()-log_bin_mean_freq.size()/5;
  int rollover_index = 0;
  int rollover_index_r2 = 0;
  vector<float> beta_vec(log_bin_mean_freq.size(),NoDataValue);
  vector<float> R_sq_vec(log_bin_mean_freq.size(),NoDataValue);
  cout << "\t\t\t\t\t 3) Now I am working on the cutoff....";
  for(int i = 0; i<cutoff;++i)
  {
    float R_sq_i,beta_i,intercept;
    vector<float> PSD_iter = slice_vector(log_bin_mean_PSD.begin()+i, log_bin_mean_PSD.end());
    vector<float> Freq_iter = slice_vector(log_bin_mean_freq.begin()+i, log_bin_mean_freq.end());
    
    least_squares_linear_regression(Freq_iter,PSD_iter,intercept,beta_i, R_sq_i);
    beta_i=-beta_i;
    beta_vec[i]=beta_i;
    R_sq_vec[i]=R_sq_i;
    cout << "log Freq: " << Freq_iter[0] << " beta: " << beta_i << " R_sq: " << R_sq_i << " " << PSD_iter.size() << "/" << log_bin_mean_PSD.size() << endl; 
    if(beta_i>beta_max)
    {
      beta_max = beta_i;
      rollover_freq = bin_mean_freq[bin_index[i]];
      rollover_index = i;
    }
    if(R_sq_i > r2_max)
    {
      r2_max = R_sq_i;
      //rollover_freq_r2 = bin_mean_freq[bin_index[i]];
      rollover_index_r2 = i;
    }
  }
  cout << "done."<< endl;
  
  cout << "Rollover index is: " << rollover_index << endl;
  cout << "rollover frequency is: " << rollover_freq << endl;
  cout << "beta_max is: " << beta_max << endl;
  
  // check if R2 gives a peak
  int n_data = int(log_bin_mean_PSD.size());
  if (rollover_index == 0 || rollover_index == n_data-1)
  {
    rollover_index = rollover_index_r2;
    cout << "I don't like that index since it will crash the program. " << endl
         << "Using R^2 to fit the rollover instead." << endl;
    cout << "The max r2 is: " << r2_max << " and index: " << rollover_index << endl; 
  }
  
  // if the rollover index is still 0, we use the entire power spectrum to fit
  if (rollover_index == 0)
  {
    rollover_index = n_data-1;
  }
  
  
  // Now get the scaling factor for the spectrum with frequencies below the 
  // rollover frequency
  
  cout << "\t\t\t\t\t 4) Slicing the vectors....";
  vector<float> PSD_below_rollover = slice_vector(log_bin_mean_PSD.begin()+2, log_bin_mean_PSD.begin()+rollover_index+1);
  vector<float> Freq_below_rollover = slice_vector(log_bin_mean_freq.begin()+2, log_bin_mean_freq.begin()+rollover_index+1);
  cout << "done."<< endl;
  
  
  
  cout << "\t\t\t\t\t 5) Doing least squares regression....";
  float beta_2,R_sq2, intercept2;
  least_squares_linear_regression(Freq_below_rollover,PSD_below_rollover,intercept2,beta_2, R_sq2);
  cout << "done."<< endl;
  
    
  // Loop through the R-squared and beta values, finding maximum in beta, and
  // point at which dR_sq/df approaches 0 to pick the rollover frequency.
//   for(int i = 0; i < beta.size();; ++i)
//   {
//     
//   }
  rollover_beta = beta_max;
  rollover_frequency = rollover_freq;
  beta = beta_vec;
  R_sq = R_sq_vec;
  sub_rollover_beta = -beta_2;
}

// calculate_background_spectrum
void LSDRasterSpectral::calculate_background_spectrum(float rollover_frequency, float beta, float log_bin_width, int N_iterations, int window_option)
{
  cout << "LINE 2224. I am going to calculate rollover frequency, so I am ignoring the one you gave to me." << endl;
  int N = RadialPSD.size();
  vector<float> background_PSD(N,0.0);
  vector<float> normalised_PSD(N,0.0);
  vector<float> ConfInt95(N,0.0);
  vector<float> normalised_ConfInt95(N,0.0);
  vector<float> normalised_ConfInt99(N,0.0);
  
  // first find rollover frequency
  float sub_rollover_beta;
  find_rollover_frequency(rollover_frequency, beta, sub_rollover_beta,log_bin_width);
  cout << "\t rollover frequency = " << rollover_frequency << "; beta = " << beta << "; sub_rollover_beta = " << sub_rollover_beta << endl;
  // Create spectra for 1000x random fractal surfaces with the same variance and
  // gradient in log-log frequency space.
  
  //--- Declarations ---
  Array2D<float> background_PSD2D(Ly,Lx,0.0);
  float chi_squared_95 = 5.991;
  float chi_squared_99 = 9.210;
  float range = 1.0;
  //int feature_order = 7;
  float RadialFreq;
  Array2D<float> SpectrumReal(Ly,Lx,0.0);
  Array2D<float> SpectrumImaginary(Ly,Lx,0.0);
  Array2D<float> fractal_raster_array(Ly,Lx,0.0);      
  Array2D<float> fractal_raster_clip(NRows,NCols,0.0);         
  Array2D<float> fractal_raster_dt(NRows,NCols,0.0);          
  Array2D<float> raster_dt(NRows,NCols,0.0);         
  Array2D<float> fractal_raster_window(NRows,NCols,0.0);
  Array2D<float> trend_plane(NRows,NCols,0.0);
  Array2D<float> window(NRows,NCols,0.0);     
  float mean_fractal,variance_fractal,mean_dt,variance_dt;
  //float variance_topo;
  //float mean_topo;
  int transform_direction;
  //mean_topo = get_mean_ignore_ndv(RasterData, NoDataValue);
  //variance_topo = get_variance_ignore_ndv(RasterData, NoDataValue, mean_topo);
  detrend2D(RasterData, raster_dt, trend_plane);
  mean_dt = get_mean_ignore_ndv(raster_dt, NoDataValue);
  variance_dt = get_variance_ignore_ndv(raster_dt, NoDataValue, mean_dt);
  //---
  cout << "\t Estimating background spectrum using fractal surfaces..." << endl;
  for(int i_iter = 0; i_iter<N_iterations; ++i_iter)
  {
    cout << flush << "\t\t\t\t\t\t\t Iteration " << i_iter+1 << " of " << N_iterations << "\r";
    // first generate a fractal surface with the required frequency scaling
    LSDRasterSpectral fractal_raster(Ly, Lx, XMinimum, YMinimum, DataResolution, NoDataValue, fractal_raster_array.copy());
    fractal_raster.rewrite_with_random_values(range);
    fractal_raster_array = fractal_raster.get_RasterData();
    
    // perform the forward fourier annalysis
    transform_direction = -1;   // this means it will be a forward transform
    fractal_raster.dfftw2D_fwd(fractal_raster_array, SpectrumReal, SpectrumImaginary,transform_direction);
    
    // Now scale spectrum
    fractal_raster.shift_spectrum(SpectrumReal, SpectrumImaginary);
    fractal_raster.scale_spectrum(SpectrumReal, SpectrumImaginary,  sub_rollover_beta/2); 
    fractal_raster.shift_spectrum_inv(SpectrumReal, SpectrumImaginary);
    
    //now perform inverse transform
    transform_direction = 1;   // this means it will be a reverse transform
    fractal_raster.dfftw2D_inv(SpectrumReal, SpectrumImaginary, fractal_raster_array, transform_direction);
    
    // Now clip the raster to the same size as the reference topography
    for (int i=0;i<NRows;++i)
    {
      for (int j=0;j<NCols;++j)
      {
        fractal_raster_clip[i][j] = fractal_raster_array[i][j];
      }
    }    
    
    // get mean and variance of the dataset
    mean_fractal = get_mean_ignore_ndv(fractal_raster_clip, NoDataValue);
    variance_fractal = get_variance_ignore_ndv(fractal_raster_clip, NoDataValue, mean_fractal);
    
    // scale dataset so that variance equals variance of real topography
    for (int i=0;i<NRows;++i)
    {
      for (int j=0;j<NCols;++j)
      {
        fractal_raster_clip[i][j] = fractal_raster_clip[i][j]*sqrt(variance_dt/variance_fractal);
      }
    }  
    mean_fractal = get_mean_ignore_ndv(fractal_raster_clip, NoDataValue);
    variance_fractal = get_variance_ignore_ndv(fractal_raster_clip, NoDataValue, mean_fractal);
    
    // rerun to remove rounding errors
    for (int i=0;i<NRows;++i)
    {
      for (int j=0;j<NCols;++j)
      {
        fractal_raster_clip[i][j] = fractal_raster_clip[i][j]*sqrt(variance_dt/variance_fractal);
      }
    }     
    mean_fractal = get_mean_ignore_ndv(fractal_raster_clip, NoDataValue);
    variance_fractal = get_variance_ignore_ndv(fractal_raster_clip, NoDataValue, mean_fractal);
    LSDRasterSpectral fractal_clip(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, fractal_raster_clip.copy());
//     fractal_clip.write_raster("test","flt");
    
    // Now get prpare for forward transform
    fractal_clip.detrend2D(fractal_raster_clip, fractal_raster_dt, trend_plane);
    fractal_clip.window_data(fractal_raster_dt, fractal_raster_window, window, window_option);
    Array2D<float> fractal_raster_padded(Ly,Lx,0.0);
    for (int i=0;i<Ly;++i)
    {
      for (int j=0;j<Lx;++j)
      {
        if (i<NRows && j<NCols)
        {
          fractal_raster_padded[i][j] = fractal_raster_window[i][j];
        }
      }
    }
    
    // perform the forward fourier analysis
    transform_direction = -1;   // this means it will be a forward transform
    fractal_clip.dfftw2D_fwd(fractal_raster_padded, SpectrumReal, SpectrumImaginary,transform_direction);
    // Now compute the periodogram of this spectrum
    fractal_clip.shift_spectrum(SpectrumReal, SpectrumImaginary);

    // CALCULATE THE DFT PERIODOGRAM
    // Multiply output by complex conjugate and normalise.
    // Note that for complex number z=x+iy, z*=x-iy, z.z* = x^2 + y^2
    fractal_clip.calculate_2D_PSD(SpectrumReal, SpectrumImaginary); 
    //float variance_ratio = 0;
    float variance_fractal_spectrum=0;
    float variance_topo_spectrum=0;
    for(int i = 0; i < P_DFT.dim1();++i)
    {
      for(int j = 0; j < P_DFT.dim2();++j)
      {
        variance_fractal_spectrum+=fractal_clip.P_DFT[i][j];
        variance_topo_spectrum+=P_DFT[i][j];
      }
    }
    // Update the average 2D spectrum
    for(int i = 0; i < Ly; ++i)
    {
      for(int j = 0; j < Lx; ++j)
      {
        background_PSD2D[i][j]+=fractal_clip.P_DFT[i][j]/float(N_iterations);
      }
    }
  }
  // Convert to radial.
  vector<float> radial_frequency(Ly*(Lx/2+1),0.0); // This is the distance from the origin in Frequency space.  Note that half of spectrum is redundant, since the fourier transform of a real dataset is symmetric, with a degeneracy of two.
  vector<float> radial_PSD(Ly*(Lx/2+1),0.0);
  int count = 0;
  for (int i=0; i < Ly; ++i)
  {
    for (int j=0; j < (Lx/2+1); ++j)
    {
      float x = j;
      float y = i;
      RadialFreq = sqrt((y - (Ly/2))*(y - (Ly/2))*dfy*dfy + (x - (Lx/2))*(x - (Lx/2))*dfx*dfx); // distance from centre to this point. Converting position in frequency into an absolute frequency
      if (RadialFreq <= NyquistFreq)  // Ignore radial frequencies greater than the Nyquist frequency as these are aliased
      {
        radial_frequency[count] = RadialFreq;
        background_PSD[count] = 2*background_PSD2D[i][j];   // Due to degeneracy
        ++count;
      }
    }
  }
    
  // Sort radial frequency
  vector<size_t> index_map;
  matlab_float_sort(radial_frequency,radial_frequency,index_map);
  // Reorder amplitudes to match sorted frequencies
  matlab_float_reorder(background_PSD,index_map,background_PSD);
  // Divide the radialPSD by the background PSD
  for(int i = 0; i<N; ++i)
  {
    normalised_PSD[i]=RadialPSD[i]/background_PSD[i];
    // 95% confidence interval using chi-squared distribution
    ConfInt95[i] = 0.5*chi_squared_95*background_PSD[i]; 
    normalised_ConfInt95[i] = 0.5*chi_squared_95;    
    normalised_ConfInt99[i] = 0.5*chi_squared_99;   
  }  
  CI95 = ConfInt95;
  normCI95 = normalised_ConfInt95;
  normCI99 = normalised_ConfInt99;
  NormalisedPSD=normalised_PSD;
  BackgroundPSD=background_PSD;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// WRAPPER FUNCTION TO CARRY OUT SPECTRAL ANALYSIS 
void LSDRasterSpectral::full_spectral_analysis(float log_bin_width, int N_iterations, int window_option)
{
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DETREND DATA
  // FIT PLANE BY LEAST SQUARES REGRESSION AND USE COEFFICIENTS TO DETERMINE
  // LOCAL SLOPE ax + by + c = z
  Array2D<float> zeta_detrend(NRows,NCols,0.0);
  Array2D<float> trend_plane(NRows,NCols,0.0);
  cout << "Detrending the data...";
  detrend2D(RasterData, zeta_detrend, trend_plane);
  cout << "finished!" << endl;

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // USE ELLIPTICAL 2D HANN (raised cosine) WINDOW ON ZETA MATRIX.
  // RETURN WINDOWED DATA AND ALSO THE SUMMED SQUARE OF THE WEIGHTING
  // COEFFICIENTS.
  Array2D<float> window(NRows,NCols,0.0);
  Array2D<float> zeta_window(NRows,NCols,0.0);
  cout << "Windowing the data...";
  window_data(zeta_detrend, zeta_window, window,window_option);
  cout << "finished!"<< endl;

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // 2D DISCRETE FAST FOURIER TRANSFORM
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // PAD DATA WITH ZEROS TO A POWER OF TWO (facilitates FFT)
  Array2D<float> zeta_padded(Ly,Lx,0.0);
  for (int i=0;i<Ly;++i)
  {
    for (int j=0;j<Lx;++j)
    {
      if (i<NRows && j<NCols)
      {
        zeta_padded[i][j] = zeta_window[i][j];
      }
    }
  }
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // DO 2D FORWARD FAST FOURIER TRANSFORM
  int transform_direction = -1;
  Array2D<float> SpectrumReal(Ly,Lx,0.0);
  Array2D<float> SpectrumImaginary(Ly,Lx,0.0);
  cout << "\t Forward transform" << endl;
  dfftw2D_fwd(zeta_padded, SpectrumReal, SpectrumImaginary, transform_direction);
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // REARRANGE SPECTRUM SO THAT ORIGIN IS AT THE CENTRE
  cout << "\t Shifting spectrum..." << endl;
  shift_spectrum(SpectrumReal, SpectrumImaginary);                                               
  
  // CALCULATE THE DFT PERIODOGRAM
  // Multiply output by complex conjugate and normalise.
  // Note that for complex number z=x+iy, z*=x-iy, z.z* = x^2 + y^2
  cout << "\t Calculating 2D periodogram..." << endl;
  float mean = get_mean_ignore_ndv(zeta_window, NoDataValue);
  float variance_hann =  get_variance_ignore_ndv(zeta_window, NoDataValue, mean);
  
  calculate_2D_PSD(SpectrumReal, SpectrumImaginary);
  // now get the total of the P_DFT
  float total_power = 0;
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if ( P_DFT[row][col] != NoDataValue)
      {
        total_power+=P_DFT[row][col];  // This calculates the sum squared amplitude, i.e. the total power, which should eqaul the total variance folllowing Parceval's theorem.
      }
    }
  }
  
  mean = get_mean_ignore_ndv(zeta_detrend, NoDataValue);
  float variance_detrend =  get_variance_ignore_ndv(zeta_detrend, NoDataValue, mean);
  mean = get_mean_ignore_ndv(zeta_padded, NoDataValue);
  float variance_padded =  get_variance_ignore_ndv(zeta_padded, NoDataValue, mean);
  cout << "variance detrend = " << variance_detrend << "; window = " << variance_hann << "; padded " << variance_padded << "; spectral power = " << total_power << "; variance_window/WSS = " << variance_hann*NRows*NCols/WSS << endl << endl << endl;
  
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // GET RADIAL POWER SPECTRUM
  // For forward transform, return the spectral power of the topography both
  // in a 2D array, and also as a one dimensional array of radial frequency
  cout << "\t Converting to radial PSD..." << endl;
  calculate_radial_PSD();

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // CALCULATE BACKGROUND FREQUENCY AND NORMALISED SPECTRUM
  cout << "\t calculate background frequency..."  << endl;
  // Initiate output vectors
  float rollover_beta = 0.1;      // This isn't used. Now the code finds a rollover beta. 
  float rollover_frequency = 0.1;     // This isn't used. Now the code finds a rollover frequency. 
  calculate_background_spectrum(rollover_frequency, rollover_beta, log_bin_width, N_iterations,window_option);
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// SPECTRUM OUTPUT
// Prints spectra, both true and normalised against background, both raw data and
// log-binned data (two separate files). 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// prints a radially averaged power spectrum, as both raw spectrum and as binned
// spectrum
void LSDRasterSpectral::print_radial_spectrum(string file_id)
{
    // First check that all the output vectors have values. If not, then fill with NoData.  Radial PSD and RadialFrequency are assumed to have been calculated.
  int N_points = RadialFrequency.size();
  vector<float> temp_points(N_points,NoDataValue);
  if(NormalisedPSD.size()==0) NormalisedPSD=temp_points;
  
  cout << "\n\t Writing radial spectrum to file" << endl;
  ofstream ofs;
  string PSD_suffix = "_radialPSD";
  string output_format = ".txt";
  string RadialPSD_file = file_id + PSD_suffix + output_format;
  ofs.open(RadialPSD_file.c_str());
  if( ofs.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << RadialPSD_file << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "Freq Wavelength PSD PSDnorm\n";
  int n_freqs = RadialFrequency.size();
  for(int i=0; i < n_freqs; ++i)
  {
    ofs << RadialFrequency[i] << " " << 1/RadialFrequency[i] << " " << RadialPSD[i] << " " << NormalisedPSD[i] << " \n";
  }
  ofs.close();
}   
// print_binned_spectrum
void LSDRasterSpectral::print_binned_spectrum(string output_id, float log_bin_width) 
{
  cout << "\n\t binning the spectrum data into logarithmically spaced bins" << endl;
    
  vector<float> bin_midpoints, bin_mean_PSD, bin_mean_freq,bin_mean_background, bin_mean_norm, bin_mean_CI, bin_CI_norm, bin_CI_norm99;
  vector< vector<float> > temp;
  log_bin_data(RadialPSD, RadialFrequency, log_bin_width, bin_midpoints, bin_mean_PSD, bin_mean_freq, temp);
  
  // First check that all the output vectors have values. If not, then fill with NoData.  Radial PSD and RadialFrequency are assumed to have been calculated.
  int N_points = RadialFrequency.size();
  int N_bins = bin_mean_freq.size();
  vector<float> temp_bins(N_bins,NoDataValue);
  vector<float> temp_points(N_points,NoDataValue);
  if(BackgroundPSD.size()==0) BackgroundPSD=temp_points;
  if(normCI95.size()==0) normCI95=temp_points;
  if(normCI99.size()==0) normCI99=temp_points;
  if(beta.size()==0) beta=temp_bins;
  if(R_sq.size()==0) R_sq=temp_bins;
  
  log_bin_data(BackgroundPSD, RadialFrequency, log_bin_width, bin_midpoints, bin_mean_background, bin_mean_freq, temp);
  log_bin_data(normCI95, RadialFrequency, log_bin_width, bin_midpoints, bin_CI_norm, bin_mean_freq, temp);
  log_bin_data(normCI99, RadialFrequency, log_bin_width, bin_midpoints, bin_CI_norm99, bin_mean_freq, temp);
  vector<int> bin_index;
  for(int i = 0; i < int(bin_mean_freq.size()); ++i)
  {
    if((bin_mean_PSD[i]!=NoDataValue)&&(bin_mean_freq[i]!=NoDataValue))
    {
      bin_index.push_back(i);
    }
  }  
  
  cout << "\t\t\t\t\t Writing log-binned spectrum file" << endl;
  // print test file
  ofstream ofs;
  // Now write this to a new file.
  string binnedPSD_suffix = "_binned_radialPSD.txt";
  string file_id = output_id;
  string binnedRadialPSD_file = file_id + binnedPSD_suffix;

  ofs.open(binnedRadialPSD_file.c_str());
  if( ofs.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << binnedRadialPSD_file << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "Freq Wavelength PSD R_sq Beta BackgroundPSD CI95 CI99 \n";
  for(int i=0; i < int(bin_index.size()); ++i)
  {
    ofs << bin_mean_freq[bin_index[i]] << " " << 1/bin_mean_freq[bin_index[i]] << " " << bin_mean_PSD[bin_index[i]] << " " << R_sq[i] << " " << beta[i] << " " << bin_mean_background[bin_index[i]] << " " << bin_CI_norm[bin_index[i]] << " " << bin_CI_norm99[bin_index[i]] << " \n";
  }
  ofs.close();
}  

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Isolate channelised portions of the landscape using an adapted method to that proposed
// by Pelletier, J. D. (2013), A robust, two-parameter method for the extraction of
// drainage networks from high-resolution digital elevation models (DEMs): Evaluation
// using synthetic and real-world DEMs, Water Resour. Res., 49, doi:10.1029/2012WR012452.
//
// This function (i) filters the DEM using a Wiener filter; (ii) calculates the
// curvature; (iii) uses quantile quantile analysis to define a curvature threshold
// for channel initiation.
//
// DTM 21/04/2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- =-=-=-=-=
LSDIndexRaster LSDRasterSpectral::IsolateChannelsWienerQQ(float area_threshold, float window_radius, string q_q_filename)
{
  cout << "\t Isolation of channelised pixels using curvature" << endl;   
  // filter
  cout << "\t\t Wiener filter" << endl;
  //float slope_percentile = 90;
  //float dt = 0.1;
  LSDRaster FilteredTopo = fftw2D_wiener();
  // calculate curvature
  vector<LSDRaster> output_rasters;
//  float window_radius = 1;
  vector<int> raster_selection(8,0.0);
  raster_selection[6]=1;
  cout << "\t\t Calculate tangential curvature" << endl;
  output_rasters = FilteredTopo.calculate_polyfit_surface_metrics(window_radius, raster_selection);
  LSDRaster curvature = output_rasters[6];
  // use q-q plot to isolate the channels
  cout << "\t\t Finding threshold using q-q plot" << endl;
  //int half_width = 100;
  LSDIndexRaster channels_init = curvature.IsolateChannelsQuantileQuantile(q_q_filename);
  // Calculate D_inf
  cout << "\t\t D_inf flow routing" << endl;
  LSDRaster Area = FilteredTopo.D_inf();
  // Reclassification of channel mask based on imposed threshold
  Array2D<int> binary_array = channels_init.get_RasterData();
  for(int i = 0; i<NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      if(channels_init.get_data_element(i,j)!=channels_init.get_NoDataValue())
      {
        if(Area.get_data_element(i,j)<=area_threshold) binary_array[i][j] = 0;
      }
    }
  }
  
  LSDIndexRaster channels(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,binary_array,GeoReferencingStrings);
  return channels;
}
LSDIndexRaster LSDRasterSpectral::IsolateChannelsWienerQQAdaptive(float area_threshold, float window_radius, string q_q_filename)
{
  cout << "\t Isolation of channelised pixels using curvature" << endl;   
  // filter
  cout << "\t\t Wiener filter" << endl;
  //float slope_percentile = 90;
  //float dt = 0.1;
  LSDRaster FilteredTopo = fftw2D_wiener();
  // calculate curvature
  vector<LSDRaster> output_rasters;
//  float window_radius = 1;
  vector<int> raster_selection(8,0.0);
  raster_selection[6]=1;
  cout << "\t\t Calculate tangential curvature" << endl;
  output_rasters = FilteredTopo.calculate_polyfit_surface_metrics(window_radius, raster_selection);
  LSDRaster curvature = output_rasters[6];
  // use q-q plot to isolate the channels
  cout << "\t\t Finding threshold using q-q plot" << endl;
  int half_width = 100;
  LSDIndexRaster channels_init = curvature.IsolateChannelsQuantileQuantileAdaptive(half_width);
  // Calculate D_inf
  cout << "\t\t D_inf flow routing" << endl;
  LSDRaster Area = FilteredTopo.D_inf();
  // Reclassification of channel mask based on imposed threshold
  Array2D<int> binary_array = channels_init.get_RasterData();
  for(int i = 0; i<NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      if(channels_init.get_data_element(i,j)!=channels_init.get_NoDataValue())
      {
        if(Area.get_data_element(i,j)<=area_threshold) binary_array[i][j] = 0;
      }
    }
  }
  
  LSDIndexRaster channels(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,binary_array,GeoReferencingStrings);
  return channels;
}


#endif
