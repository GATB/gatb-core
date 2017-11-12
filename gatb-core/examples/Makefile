# example GATB-core makefile

.PHONY: all clean

# change this to the folder where gatb-core is
GATB=$(shell pwd)/../

GATB_LIB=$(GATB)/build/lib

CXXFLAGS =  -std=c++0x -O3  -I$(GATB)/src  -I$(GATB)/build/include -I$(GATB)/thirdparty 
LDFLAGS=   -L$(GATB_LIB) -lgatbcore  -lpthread -lz -lhdf5 -std=c++0x -ldl -static


SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%,$(SRCS))

all: $(OBJS)

clean:
	rm -fv $(OBJS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LDFLAGS)

