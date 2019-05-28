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
	LSDRaster DEM("/data/Geog-c2s2/toku/final", DEM_Format);

  //Fill
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions, FilledDEM);

  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  //get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, 10000);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

  vector<int> basin_junctions = ChanNetwork.ExtractBasinJunctionOrder(7, FlowInfo);

  for (int i = 0; i < int(basin_junctions.size()); ++i){

    vector<int> sub_basin_sources = ChanNetwork.get_all_source_nodes_of_an_outlet_junction(basin_junctions[i]);

    LSDJunctionNetwork SubChanNetwork(sub_basin_sources, FlowInfo);

    LSDIndexRaster SOArray = SubChanNetwork.StreamOrderArray_to_LSDIndexRaster();

    stringstream ss;
    ss << "/data/Geog-c2s2/toku/sub_net_" << i;

    SOArray.write_raster(ss.str(), DEM_Format);

    LSDStrahlerLinks Links(SubChanNetwork, FlowInfo);
    Links.CalculateTokunagaIndexes(SubChanNetwork, FlowInfo);

    stringstream ss2;
    ss2 << "toku_subnet_" << i;
    Links.WriteTokunagaData("/data/Geog-c2s2/toku/", ss2.str());


    stringstream ss3;
    ss3 << "/data/Geog-c2s2/toku/toku_orders_" << i;
    LSDIndexRaster source_data = Links.WriteTokunagaRaster(FlowInfo);
    source_data.write_raster(ss3.str(), DEM_Format);

  }

  LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  SOArray.write_raster("/data/Geog-c2s2/toku/main_net", DEM_Format);

  float hs_azimuth = 315;
  float hs_altitude = 45;
  float hs_z_factor = 1;
  LSDRaster hs_raster = FilledDEM.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

  FilledDEM.write_raster("/data/Geog-c2s2/toku/hillshade", DEM_Format);

}
