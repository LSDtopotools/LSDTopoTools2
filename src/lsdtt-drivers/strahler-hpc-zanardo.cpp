
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

  // stringstream ss3;
  // ss3 << Outpath << "full_network_" << DEMname;
  // ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, ss3.str());

  int nearest_node = atoi(argv[4]);
  std::cout << "\nNearest node: " << nearest_node << '\n';

  int nearest_junction = ChanNetwork.get_Junction_of_Node(nearest_node, FlowInfo);
  std::cout << "\nConverted to junction:" << nearest_junction << '\n';

  int junct_node = ChanNetwork.get_Node_of_Junction(nearest_junction);

  std::cout << "\nback to node: " << junct_node << '\n';

  int blah = ChanNetwork.get_Junction_of_Node(junct_node, FlowInfo);

  std::cout << "\nback to junction: " << blah << '\n';

  vector<int> sub_basin_sources = ChanNetwork.get_all_source_nodes_of_an_outlet_junction(nearest_junction);

  std::cout << "Source count: " << sub_basin_sources.size() << '\n';

  LSDJunctionNetwork SubChanNetwork(sub_basin_sources, FlowInfo);

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
