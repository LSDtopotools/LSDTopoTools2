CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -Wall -O3
SOURCES = zanardo-2.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDShapeTools.cpp \
    ../LSDStatsTools.cpp \
		../LSDJunctionNetwork.cpp \
		../LSDIndexChannel.cpp \
		../LSDChannel.cpp \
		../LSDMostLikelyPartitionsFinder.cpp \
	  ../LSDStrahlerLinks.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=zanardo-2.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
