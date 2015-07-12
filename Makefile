# /********************************************************************************************
# * File:		Makefile
# * Author:		$LastChangedBy: bolitho $
# * Revision:		$Revision: 82 $
# * Last Updated:	$LastChangedDate: 2005-09-26 07:31:28 -0400 (Mon, 26 Sep 2005) $
# ********************************************************************************************/

TARGET=PoissonRecon
SOURCE=CmdLineParser.cpp Factor.cpp Geometry.cpp MarchingCubes.cpp ply.cpp plyfile.cpp Time.cpp MultiGridOctest.cpp

CFLAGS += -fpermissive
LFLAGS +=

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -pg
LFLAGS_RELEASE = -O3 -pg

SRC = ./
BIN = Bin/Linux/
INCLUDE = /usr/include/

CC=gcc
CXX=g++

OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(SOURCE))))

all: CFLAGS += $(CFLAGS_DEBUG)
all: LFLAGS += $(LFLAGS_DEBUG)
all: $(BIN)$(TARGET)

release: CFLAGS += $(CFLAGS_RELEASE)
release: LFLAGS += $(LFLAGS_RELEASE)
release: $(BIN)$(TARGET)

clean:
	rm -f $(BIN)$(TARGET)
	rm -f $(OBJECTS)

$(BIN)$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LFLAGS)

$(BIN)%.o: $(SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

