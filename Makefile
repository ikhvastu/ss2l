CC=g++
CFLAGS=
LDFLAGS=`root-config --libs --cflags`
SOURCES=readTree.cc mt2_bisect.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=readTree

all:
	$(CC) $(CFLAGS) $(SOURCES) $(LDFLAGS) -o $(EXECUTABLE)

clean:
	rm -rf *o $(EXECUTABLE)