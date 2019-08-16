
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
#include "../LSDBasin.hpp"

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

  // Get topo derivatives
  vector<int> raster_selection;

  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(1); //aspect
  raster_selection.push_back(0); //curvature
  raster_selection.push_back(0); //plan curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  int WindowSize = 9;

  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(WindowSize, raster_selection);
  LSDRaster Slope = Surfaces[1];
  LSDRaster Aspect = Surfaces[2];

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions, FilledDEM);

  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  //get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, 103);

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

  int max_order = ChanNetwork.get_maximum_stream_order();

  if (max_order > 5){

    vector<int> basin_junctions = ChanNetwork.ExtractBasinJunctionOrderKeepEdgeBasins(max_order - 1, FlowInfo);

    std::cout << "Number of basins to process: " << basin_junctions.size() << '\n';

    for (int i = 0; i < int(basin_junctions.size()); ++i){

      vector<int> sub_basin_sources = ChanNetwork.get_all_source_nodes_of_an_outlet_junction(basin_junctions[i]);

      LSDJunctionNetwork SubChanNetwork(sub_basin_sources, FlowInfo);

      // LSDBasin sub_basin = LSDBasin(basin_junctions[i], FlowInfo, SubChanNetwork);
      //
      // sub_basin.set_AspectMean(FlowInfo, Aspect);
      // sub_basin.set_SlopeMean(FlowInfo, Slope);
      // sub_basin.set_DrainageDensity();
      // float max_elev = sub_basin.CalculateBasinMax(FlowInfo, FilledDEM);
      // float min_elev = sub_basin.CalculateBasinMin(FlowInfo, FilledDEM);
      // int n_heads = int(sub_basin_sources.size());
      //
      // stringstream ss;
      // ss << Outpath << "BasinData_" << DEMname << "_" << i << ".csv";
      //
      // ofstream basin_data_writer;
      // basin_data_writer.open(ss.str().c_str());
      // basin_data_writer << n_heads << "," << sub_basin.get_SlopeMean() << "," << sub_basin.get_AspectMean() << "," << sub_basin.get_Area() << "," << sub_basin.get_DrainageDensity() << "," << max_elev << "," << min_elev << endl;
      // basin_data_writer.close();

    }
}

else{

  std::cout << "There are no suitable basins to process, the max order found is: " << max_order << '\n';

}

}
