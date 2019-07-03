
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
	if (nNumberofArgs!=6)
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

  stringstream ss3;
  ss3 << Outpath << "full_network_" << DEMname;
  ChanNetwork.StreamOrderArray_to_WGS84CSV(ss3.str());

  double Lat = atof(argv[4]);
  double Long = atof(argv[5]);
  LSDCoordinateConverterLLandUTM converter;

  double northing, easting;
  int UTM_zone;
  bool is_North;
  FlowInfo.get_UTM_information(UTM_zone, is_North);
  int eId = 22;             // defines the ellipsiod. This is WGS
  converter.LLtoUTM_ForceZone(eId, Lat, Long, northing, easting, UTM_zone);

  std::cout << "Searching for nearest channel" << '\n';
  std::cout << FlowInfo.get_NRows() << '\n';
  std::cout << FlowInfo.get_NCols() << '\n';
  std::cout << "------------------" << '\n';
  int nearest_node = ChanNetwork.get_nodeindex_of_nearest_channel_for_specified_coordinates(easting, northing, 6, 25, FlowInfo);
  int nearest_junction = ChanNetwork.get_Junction_of_Node(nearest_node, FlowInfo);

  vector<int> sub_basin_sources = ChanNetwork.get_all_source_nodes_of_an_outlet_junction(nearest_junction);

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
