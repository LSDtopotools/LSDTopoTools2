
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
	if (nNumberofArgs!=5)
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

  int threshold = 103;  // needs to be 0.09 km2 in pixels

  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  std::cout << "Built junction network" << '\n';

  int nearest_node = atoi(argv[4]);
  std::cout << "\nNearest node: " << nearest_node << '\n';

  int nearest_junction = ChanNetwork.find_upstream_junction_from_channel_nodeindex(nearest_node, FlowInfo);
  std::cout << "\nConverted to junction:" << nearest_junction << '\n';

  vector<int> sub_basin_sources = ChanNetwork.get_all_source_nodes_of_an_outlet_junction(nearest_junction);

  LSDJunctionNetwork SubChanNetwork(sub_basin_sources, FlowInfo);

  int max_order = ChanNetwork.get_maximum_stream_order();

  if (max_order > 6){

    vector<int> basin_junctions = SubChanNetwork.ExtractBasinJunctionOrderKeepEdgeBasins(7, FlowInfo);

    std::cout << "Number of junctions: " << basin_junctions.size() << '\n';

    if (int(basin_junctions.size()) > 0){

      for (int i = 0; i < int(basin_junctions.size()); ++i){

        vector<int> sub_basin_sources = SubChanNetwork.get_all_source_nodes_of_an_outlet_junction(basin_junctions[i]);

        LSDJunctionNetwork SubSubChanNetwork(sub_basin_sources, FlowInfo);

        LSDStrahlerLinks Links(SubSubChanNetwork, FlowInfo);
        Links.CalculateTokunagaIndexes(SubSubChanNetwork, FlowInfo);
        Links.calculate_lengths(FlowInfo);

        stringstream ss;
        ss << DEMname << "_" << i;
        Links.WriteTokunagaData(Outpath, ss.str());

        stringstream ss2;
        ss2 << Outpath << "toku_network_" << DEMname << "_" << i;
        Links.WriteTokunagaChannelsCSV(SubSubChanNetwork, ss2.str());
      }
    }
    else{

      LSDStrahlerLinks Links(SubChanNetwork, FlowInfo);
      Links.CalculateTokunagaIndexes(SubChanNetwork, FlowInfo);
      Links.calculate_lengths(FlowInfo);

      stringstream ss;
      ss << DEMname << "_" << 1;
      Links.WriteTokunagaData(Outpath, ss.str());

      stringstream ss2;
      ss2 << Outpath << "toku_network_" << DEMname << "_" << 1;
      Links.WriteTokunagaChannelsCSV(SubChanNetwork, ss2.str());

    }


  } else{

    LSDStrahlerLinks Links(SubChanNetwork, FlowInfo);
    Links.CalculateTokunagaIndexes(SubChanNetwork, FlowInfo);
    Links.calculate_lengths(FlowInfo);

    stringstream ss;
    ss << DEMname << "_" << 1;
    Links.WriteTokunagaData(Outpath, ss.str());

    stringstream ss2;
    ss2 << Outpath << "toku_network_" << DEMname << "_" << 1;
    Links.WriteTokunagaChannelsCSV(SubChanNetwork, ss2.str());

  }

}
