
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

  int threshold = 103;  // needs to be 0.09 km2 in pixels

  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  std::cout << "Built junction network" << '\n';

  stringstream ss;
  ss << Outpath << "full_network_" << DEMname;
  ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, ss.str());

}
