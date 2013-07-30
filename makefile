CC = g++
CFLAGS = -c -Wall -g
LDFLAGS = 
SOURCES = testFlow.cpp flow.cpp ObsParameters.cpp lambda.cpp 
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = testFlow

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
