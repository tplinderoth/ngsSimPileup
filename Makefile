CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall
CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)
BIN = ngsSimPileup

all: $(BIN)

.PHONY: clean

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp > $*.d

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) *.o

clean:
	rm -f $(BIN) *.o *.d
