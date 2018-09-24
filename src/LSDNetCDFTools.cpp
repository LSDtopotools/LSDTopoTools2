// LSDNetCDFTools.cpp

/// Note: you need to compile with -lnetcdf and -lnetcdf_cxx4 flags
/// You will need the need the netCDF-4 C++ libraries installed
/// Most Linux package managers will install this for you using
/// yum, dnf, or apt-get.

/// A tool for reading netCDF-4 files into LSDTopoTools
///
/// NetCDF (The network common data format is a file format for storing
/// multivariate climate, atmospheric, etc. data. (Although it can be used to store
/// any type of gridded data. You could use it to store states of topography through
/// time, for example.)

/// Full documentation of the netCDF C++ API can be found at:
/// http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-cxx
///
/// @author DAV
/// @date 2016-11-05

#include <iostream>
#include <netcdf>
#include "TNT/tnt.h"

using namespace netCDF;
using namespace netCDF::exceptions;

// Return this in event of a problem.
static const int NC_ERR = 2;

// Originally had it as a class but not much point at this stage
// just make it a namespace instead
namespace LSDNetCDFTools
{

/// @brief Reads in a netCDF field and writes it into
/// a TNT array 1D.
/// @return An error code, (0 for success, NC_ERR for failure)
/// @details Template allows to specify different TNT array data types
/// e.g. double, float, int etc.
/// @param reference to a tnt::array1d, the name of a netcdf file to read
/// @note 1d not yet implemented...see 2d.
template<typename T>
int read_tnt_array1d(TNT::Array1D<T>*const tnt_array);

/// @brief Reads in a netCDF field and writes it into
/// a TNT array 2D.
/// @return An error code, (0 for success, NC_ERR for failure)
/// @details Template allows to specify different TNT array data types
/// e.g. double, float, int etc.
/// @param reference to a tnt::array2d, the name of a netcdf file to read
template<typename T>
int read_tnt_array2d(TNT::Array2D<T>*const tnt_array, std::string nc_file);
// The pointer to the tnt array is constant, but not the actual array itself
// i.e. we can change the array thtough the pointer, but ensure the pointer
// remains constant and cannot be changed to point to something else, which would be 
// a bad idea...

/// @brief Reads in a netCDF field and writes it into
/// a TNT array 3D.
/// @return An error code, (0 for success, NC_ERR for failure)
/// @details Template allows to specify different TNT array data types
/// e.g. double, float, int etc.
/// @param reference to a tnt::array3d, the name of a netcdf file to read
template<typename T>
int read_tnt_array3d(TNT::Array3D<T>*const tnt_array, std::string nc_file);

//=-=-=-=-=-=-=-=-=
// IMPLEMENTATIONS
//=-=-=-=-=-=-=-=-=

template<typename T>
int read_tnt_array3d(TNT::Array3D<T>*const tnt_array, std::string nc_file)
{
  int n_timesteps = tnt_array->dim1();
  int rows = tnt_array->dim2();
  int cols = tnt_array->dim3();

  std::cout << "Timesteps: " << n_timesteps << ", " << "Rows: " << rows <<
               ", " << "Cols: " << cols << std::endl;

  // Should introduce some error checking here to make sure sure
  // the dims in the netcdf file match the referenced array passed...
  try
  {
  // Open the file.
  NcFile dataFile(nc_file, NcFile::read);

  // Get the latitude and longitude variables and read data.
  NcVar latVar, lonVar;

//  latVar = dataFile.getVar("latitude");
//  if(latVar.isNull()) return NC_ERR;
//  lonVar = dataFile.getVar("longitude");
//  if(lonVar.isNull()) return NC_ERR;
//  lonVar.getVar(lons);
//  latVar.getVar(lats);

  NcVar presVar;
  presVar = dataFile.getVar("pressure");
  if(presVar.isNull()) return NC_ERR;

  // Read in data to array
  // Works because the elements are contiguous
  // in memory.
  presVar.getVar(tnt_array[0][0][0]);

  for (size_t timestep = 0; timestep < n_timesteps; timestep++)
  {
  //for (int timestep = 0; timestep < NREC; timestep++)
    std::cout << "TIMESTEP: " << timestep<< std::endl;
    for (int row = 0; row < rows; row++)
    {
      for (int col = 0; col < cols; col++)
      {
        // derefernce our pointer to array with the (*ptr)[][] notation...
        std::cout << (*tnt_array)[timestep][row][col] << " ";
        //if((*tnt_array)[rec][lat][lon] != (float) (SAMPLE_PRESSURE + i)) return NC_ERR;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

  } // next record

  // The file is automatically closed by the destructor. This frees
  // up any internal netCDF resources associated with the file, and
  // flushes any buffers.

  std::cout << "*** SUCCESS reading example file " << nc_file << std::endl;
  return 0;

  }
  catch(NcException& e)
  {
     e.what();
     std::cout<<"FAILURE**************************"<<std::endl;
     return NC_ERR;
  }
}

template<typename T>
int read_tnt_array1d(TNT::Array1D<T>*const tnt_array)
{
  //1d should be easy
}

template<typename T>
int read_tnt_array2d(TNT::Array2D<T>*const tnt_array, std::string nc_file)
{
   int NX = tnt_array->dim1();
   int NY = tnt_array->dim2();

   std::cout << "Rows: " << NX << ", " << "Cols: " << NY << std::endl;

   try
   {
   // This is a C-style array.
   // We tend not to use these in LSDTopoTools
   // but typically used in netCDF.
   // int dataIn[NX][NY];

   // This is a TNT (Template Numerical Toolkit)
   // 2D array (allocated at runtime)
   //tnt_array(NX,NY,0);

   // Open the file for read access
   NcFile dataFile(nc_file, NcFile::read);

   // Retrieve the variable named "data"
   NcVar data = dataFile.getVar("data");
   if(data.isNull()) return NC_ERR;

   // Ingest the data into the array
   // We pass a reference to the first element of
   // the TNT array [0][0]. getVar writes in the values.
   data.getVar(tnt_array[0][0]);

   // Note: with a C-style array you would do this:
   // data.getVar(dataIn)
   // Since C arrays are already pointers to memory.

   // Check the values
   for (int i = 0; i < NX; i++)
   {
      for (int j = 0; j < NY; j++)
      {
         // Check that we haven't gone out of bounds...
         //if (*tnt_array[i][j] != i * NY + j) return NC_ERR;
         //else
         {
           std::cout << (*tnt_array)[i][j] << " ";
         }
      }
      std::cout << std::endl;
   }
   
   // The netCDF file is automatically closed by the NcFile destructor
   std::cout << "*** SUCCESS reading example file " << nc_file << std::endl;

   return 0;
   }
   catch(NcException& e)
   {
     e.what();
     std::cout<<"FAILURE*************************************"<< std::endl;
     return NC_ERR;
   }
}
}

// Test.

//using namespace LSDNetCDFTools;

int main()
{
  // Test reading in a 2D netcdf file
  TNT::Array2D<int> myArray2D(6,12, 0);

  LSDNetCDFTools::read_tnt_array2d(&myArray2D, "simple_xy.nc");

  // Test reading in a 3D (i.e. one var + timestep)
  TNT::Array3D<float> myArray3D(2, 6,12, 0.0);

  LSDNetCDFTools::read_tnt_array3d(&myArray3D, "pres_temp_3D.nc");
}
