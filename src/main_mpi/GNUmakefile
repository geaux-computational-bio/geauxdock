include ../GNUmakefile.common.mk


EXE = server client_dummy client_cpu client_mic client_gpu
OBJ += load.o data.o rmsd.o util.o util_print.o


default: $(EXE)


server: server.o $(OBJ) util_optimize_soa_cpu.o
	$(CXX_MPI) $(LIBPATH) $^ -o $@ $(LINKFLAGS)




client_dummy: client_dummy.o
	$(CXX_MPI) $(LIBPATH) $^ -o $@ $(LINKFLAGS)

client_cpu: client_cpu.o backend_cpu.o $(OBJ)
	$(CXX_MPI) $(LIBPATH) $^ -o $@ $(LINKFLAGS)

client_mic: client_mic.o backend_mic.o $(OBJ)
	$(CXX_MPI) $(LIBPATH) $^ -o $@ $(LINKFLAGS)

client_gpu: client_gpu.o prepair.o backend_gpu.o $(OBJ)
	$(CXX_MPI) $(LIBPATH) $^ -o $@ $(LINKFLAGS) -lcudart




client_cpu.o: client_cpu_mic.C
	$(CXX_MPI) $(CXXFLAGS_HOST) -c $< -o $@

client_mic.o: client_cpu_mic.C
	$(CXX_MPI) $(CXXFLAGS_HOST) -c $< -o $@

client_gpu.o: client_gpu.C
	$(CXX_MPI) $(CXXFLAGS_HOST) $(MARCRO_GPU) -c $< -o $@

%.o: %.C
	$(CXX_MPI) $(CXXFLAGS_HOST) -c $< -o $@



include ../GNUmakefile.common.rules

