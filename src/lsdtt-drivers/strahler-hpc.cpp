
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

  //Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number of inputs.";
		exit(EXIT_FAILURE);
	}

  //get input args
  string Inpath = argv[1];
  string Outpath = argv[2];
  string DEMname = argv[3];

  string DEM_Format = "bil";

  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

	//load the DEM
	LSDRaster DEM(Inpath + DEMname, DEM_Format);

  //Fill
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions, FilledDEM);

  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  //get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, 3500);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

  int max_order = ChanNetwork.get_maximum_stream_order();

  if (max_order > 5){

    vector<int> basin_junctions = ChanNetwork.ExtractBasinJunctionOrderKeepEdgeBasins(max_order - 1, FlowInfo);

    std::cout << "Number of basins to process: " << basin_junctions.size() << '\n';

    for (int i = 0; i < int(basin_junctions.size()); ++i){

      vector<int> sub_basin_sources = ChanNetwork.get_all_source_nodes_of_an_outlet_junction(basin_junctions[i]);

      LSDJunctionNetwork SubChanNetwork(sub_basin_sources, FlowInfo);

      LSDStrahlerLinks Links(SubChanNetwork, FlowInfo);
      Links.CalculateTokunagaIndexes(SubChanNetwork, FlowInfo);
      Links.calculate_lengths(Flowinfo);

      stringstream ss;
      ss << DEMname << "_" << i;
      Links.WriteTokunagaData(Outpath, ss.str());


      stringstream ss2;
      ss2 << Outpath << "toku_network_" << DEMname << "_" << i;
      Links.WriteTokunagaChannelsCSV(SubChanNetwork, ss2.str());

    }

    stringstream ss3;
    ss3 << Outpath << "full_network_" << DEMname;
    ChanNetwork.StreamOrderArray_to_WGS84CSV(ss3.str());

}

else{

  std::cout << "There are no suitable basins to process, the max order found is: " << max_order << '\n';

}

}
