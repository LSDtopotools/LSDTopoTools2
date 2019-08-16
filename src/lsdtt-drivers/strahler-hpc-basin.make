CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -Wall -O3
SOURCES = strahler-hpc-basin.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDShapeTools.cpp \
    ../LSDStatsTools.cpp \
		../LSDJunctionNetwork.cpp \
		../LSDIndexChannel.cpp \
		../LSDChannel.cpp \
		../LSDMostLikelyPartitionsFinder.cpp \
	  ../LSDStrahlerLinks.cpp \
		../LSDBasin.cpp \
		../LSDCRNParameters.cpp \
		../LSDSpatialCSVReader.cpp \
		../LSDParticle.cpp
		
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=strahler-hpc-basin.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
