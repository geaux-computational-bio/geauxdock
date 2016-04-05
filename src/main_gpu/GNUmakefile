include ../GNUmakefile.common.mk


EXE := geauxdock
OBJ += geauxdock.o
OBJ += prepair.o backend_gpu.o
OBJ += util_optimize_soa_gpu.o
OBJ += load.o data.o rmsd.o util.o util_print.o


default: $(EXE)


$(EXE): $(OBJ)
	$(CXX_HOST) $(LIBPATH) $(OBJ) -o $@ $(LINKFLAGS) -lcudart

geauxdock.o: ./geauxdock.C
	$(CXX_HOST) $(CXXFLAGS_HOST) $(MARCRO_GPU) -c $< -o $@


#NT := $(GD) # number of temperatures passing to the runtime
NT := 2 # number of temperatures passing to the runtime
include ../GNUmakefile.common.rules

