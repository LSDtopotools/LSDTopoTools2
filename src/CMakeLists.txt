cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(LSDTopoTools)

# Dependencies
# ============
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries     (PkgConfig::FFTW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof -O3 -std=c++11")

# We test if the open cv library is found. 
# If not we compile without lsdtt-valley-width program
# OpenCV
find_package(OpenCV)
if(OpenCV_FOUND)
    message("I found OpenCV and am compiling the floodplain code")
    include_directories( ${OpenCV_INCLUDE_DIRS})

    # Build
    # =====
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
    LSDParticleColumn.hpp
    LSDRaster.hpp
    LSDRasterAggregator.hpp
    LSDRasterInfo.hpp
    LSDRasterMaker.hpp
    LSDRasterSpectral.hpp
    LSDSoilHydroRaster.hpp
    LSDStatsTools.hpp
    LSDShapeTools.hpp
    LSDSpatialCSVReader.hpp
    LSDStrahlerLinks.hpp
    LSDFloodplain.hpp
    LSDTerrace.hpp
    )

    string(REGEX REPLACE "([^;]+).hpp" "\\1.cpp" LIB_SOURCES "${LIB_HEADERS}")

    set(EXE_SOURCES
    lsdtt-drivers/lsdtt-basic-metrics.cpp
    lsdtt-drivers/lsdtt-channel-extraction.cpp
    lsdtt-drivers/lsdtt-chi-mapping.cpp
    lsdtt-drivers/lsdtt-cosmo-tool.cpp
    lsdtt-drivers/lsdtt-hillslope-channel-coupling.cpp
    lsdtt-drivers/lsdtt-valley-metrics.cpp
    )

    set(EXE_FILES)

    add_library(${PROJECT_NAME} ${LIB_SOURCES} ${LIB_HEADERS} )

    foreach(src_file ${EXE_SOURCES})
    string(REGEX REPLACE "([^;]+)/([^;]+).cpp" "\\2" exe_file ${src_file})
    list(APPEND EXE_FILES ${exe_file})
    add_executable(${exe_file} ${src_file})
    target_link_libraries(${exe_file} ${PROJECT_NAME} PkgConfig::FFTW)
    target_link_libraries(${exe_file} ${PROJECT_NAME} ${OpenCV_LIBS} )
    endforeach(src_file ${EXE_SOURCES})


    # Install
    # =======
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(
    TARGETS ${EXE_FILES} ${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else()
    message("I am compiling without the openCV library")
    
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
    LSDParticleColumn.hpp
    LSDRaster.hpp
    LSDRasterAggregator.hpp
    LSDRasterInfo.hpp
    LSDRasterMaker.hpp
    LSDRasterSpectral.hpp
    LSDSoilHydroRaster.hpp
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
    TARGETS ${EXE_FILES} ${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
  endif()
  
