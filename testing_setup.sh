#!/bin/bash

cd tests/
mkdir results
cd fixtures
wget https://github.com/LSDtopotools/ExampleTopoDatasets/raw/master/coweeta.bil
wget https://github.com/LSDtopotools/ExampleTopoDatasets/raw/master/coweeta.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_ASPECT.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_ASPECT.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_CLASS.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_CLASS.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_CURV.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_CURV.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_PLFMCURV.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_PLFMCURV.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_PROFCURV.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_PROFCURV.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_SLOPE.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_SLOPE.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_SMOOTH.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_SMOOTH.hdr
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_TANCURV.bil
wget https://github.com/LSDtopotools/RegressionData/raw/master/coweeta_output_TANCURV.hdr
cd ..
