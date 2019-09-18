cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(LSDTopoTools)

# Dependencies
# ============

# TODO: issue FFTW3LibraryDepends.cmake missing in most fftw3 installs
# https://github.com/FFTW/fftw3/issues/130
#find_package(fftw3 REQUIRED)
#message(STATUS "Found fttw3: ${FFTW3_LIBRARY_DIRS}/fftw3")

# TODO: check which components are needed (+ not working with pcl conda-forge)
#find_package(pcl REQUIRED COMPONENTS)
#message(STATUS "Found pcl: ${PCL_LIBRARY_DIRS}/pcl")


# Build
# =====
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof")

set(LIB_HEADERS
  LSDBasin.hpp
  LSDChiTools.hpp
  LSDChiNetwork.hpp
  LSDChannel.hpp
  LSDCosmoData.hpp
  LSDCosmoRaster.hpp
  LSDCRNParameters.hpp
  LSDFlowInfo.hpp
  LSDIndexChannel.hpp
  LSDIndexChannelTree.hpp
  LSDIndexRaster.hpp
  LSDJunctionNetwork.hpp
  LSDMostLikelyPartitionsFinder.hpp
  LSDParameterParser.hpp
  LSDParticle.hpp
  LSDRaster.hpp
  LSDRasterAggregator.hpp
  LSDRasterInfo.hpp
  LSDRasterMaker.hpp
  LSDRasterSpectral.hpp
  LSDStatsTools.hpp
  LSDShapeTools.hpp
  LSDSpatialCSVReader.hpp
  LSDStrahlerLinks.hpp
)

string(REGEX REPLACE "([^;]+).hpp" "\\1.cpp" LIB_SOURCES "${LIB_HEADERS}")

set(EXE_SOURCES
  lsdtt-drivers/lsdtt-basic-metrics.cpp
  lsdtt-drivers/lsdtt-channel-extraction.cpp
  lsdtt-drivers/lsdtt-chi-mapping.cpp
  lsdtt-drivers/lsdtt-cosmo-tool.cpp
  lsdtt-drivers/lsdtt-hillslope-channel-coupling.cpp
  lsdtt-drivers/lsdtt-river-clusters.cpp
)

set(EXE_FILES)

add_library(${PROJECT_NAME} ${LIB_SOURCES} ${LIB_HEADERS} )

foreach(src_file ${EXE_SOURCES})
  string(REGEX REPLACE "([^;]+)/([^;]+).cpp" "\\2" exe_file ${src_file})
  list(APPEND EXE_FILES ${exe_file})
  add_executable(${exe_file} ${src_file})
  target_link_libraries(${exe_file} ${PROJECT_NAME} fftw3)
endforeach(src_file ${EXE_SOURCES})


# Install
# =======

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(
  PROGRAMS ${EXE_FILES}
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)
