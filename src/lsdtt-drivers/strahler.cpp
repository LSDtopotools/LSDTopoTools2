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
	LSDRaster DEM("/LSDTopoTools/EGU-Workshop/data/Pozo_1m/Pozo_DTM", DEM_Format);

  //Fill
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions, FilledDEM);

  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  //get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, 2000);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

  vector<int> basin_junctions = ChanNetwork.ExtractBasinJunctionOrderKeepEdgeBasins(4, FlowInfo);

  for (int i = 0; i < int(basin_junctions.size()); ++i){

    vector<int> sub_basin_sources = ChanNetwork.get_all_source_nodes_of_an_outlet_junction(basin_junctions[i]);

    LSDJunctionNetwork SubChanNetwork(sub_basin_sources, FlowInfo);

    LSDIndexRaster SOArray = SubChanNetwork.StreamOrderArray_to_LSDIndexRaster();

    stringstream ss;
    ss << "/LSDTopoTools/strahler/sub_net_" << i;

    SOArray.write_raster(ss.str(), DEM_Format);

    LSDStrahlerLinks Links(SubChanNetwork, FlowInfo);
    Links.CalculateTokunagaIndexes(SubChanNetwork, FlowInfo);

    stringstream ss2;
    ss2 << "toku_subnet_" << i;
    Links.WriteTokunagaData("/LSDTopoTools/strahler/", ss2.str());


    stringstream ss3;
    ss3 << "/LSDTopoTools/strahler/toku_orders_" << i;
    LSDIndexRaster source_data = Links.WriteTokunagaRaster(FlowInfo);
    source_data.write_raster(ss3.str(), DEM_Format);

  }

  LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  SOArray.write_raster("/LSDTopoTools/strahler/main_net", DEM_Format);

  float hs_azimuth = 315;
  float hs_altitude = 45;
  float hs_z_factor = 1;
  LSDRaster hs_raster = FilledDEM.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

  hs_raster.write_raster("/LSDTopoTools/strahler/hillshade", DEM_Format);

}