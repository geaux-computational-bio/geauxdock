
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
EXE = $(SRC:.cpp=)

CXX = g++ -O2 -std=c++14


LINKFLAGS += -lyeahc -lrt


all: $(EXE)

% : %.cpp
	$(CXX) $(HEADPATH) $(LIBPATH) $< -o $@ $(LINKFLAGS)

echo:
#	@(echo $(HEADERS))
	@(echo $(SRC))
	@(echo $(OBJ))
	@(echo $(EXE))
clean:
	$(RM) $(OBJ)
	$(RM) $(EXE)

