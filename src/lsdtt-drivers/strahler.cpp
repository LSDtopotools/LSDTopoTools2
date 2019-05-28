//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Compile_Dreich_lengths.cpp
//
// Driver created to generate channel length measurements data at a range of input resolutions.
//
// Driver expects unfilled DEMs at a range of grid resolutions in the input directory
//
// Run the driver with the following arguments:
//
// path to the DEM files with a trailing slash
// output path with trailing slash
// DEM filename Prefix
// file extension without the dot
//
// DEM naming convention should be <prefix>_<resolution>_DEM
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Stuart W.D. Grieve
// University of Edinburgh
// November 2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDStrahlerLinks.hpp"

int main(int nNumberofArgs, char *argv[])
{

  string DEM_Format = "bil";

  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

	//load the DEM
	LSDRaster DEM("/LSDTopoTools/strahler/srtm", DEM_Format);

  //Fill
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions, FilledDEM);

  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  //get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, 200);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

  LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  SOArray.write_raster("/LSDTopoTools/strahler/srtm_so", DEM_Format);

  LSDStrahlerLinks Links(ChanNetwork, FlowInfo);

  Links.CalculateTokunagaIndexes(ChanNetwork, FlowInfo);

  Links.WriteTokunagaData("/LSDTopoTools/strahler/", "srtm");

  LSDIndexRaster source_data = Links.WriteTokunagaRaster(FlowInfo);

  source_data.write_raster("/LSDTopoTools/strahler/srtm_t_order", DEM_Format);
}
