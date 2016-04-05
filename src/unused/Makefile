# makefile constants

TARGET_GPU := 0
TARGET_CPU := 1
TARGET_MIC := 2





HOST := gpuhost

ifeq ($(HOST), cpuhost)
# module load impi/4.1.3.048/intel64
        LIBPATH += -L/usr/lib/x86_64-redhat-linux5E/lib64
	TARGET_DEVICE := $(TARGET_CPU)

else ifeq ($(HOST), michost)
	TARGET_DEVICE := $(TARGET_MIC)

else ifeq ($(HOST), gpuhost)
	LIBPATH += -L/usr/local/packages/cuda/6.5/lib64
	TARGET_DEVICE := $(TARGET_GPU)
	GPU := K20X
        NGPU := 1

else ifeq ($(HOST), shelob)
	LIBPATH += -L/usr/local/packages/cuda/5.5.22/lib64
	TARGET_DEVICE := $(TARGET_GPU)
	GPU := K20X
        NGPU := 1


else ifeq ($(HOST), ece)
	LIBPATH := -L/usr/local/cuda-6.5/lib64
	TARGET_DEVICE := $(TARGET_GPU)
	GPU := GTX780
        NGPU := 1

else ifeq ($(HOST), lasphi)
	LIBPATH := -L/usr/local/cuda-7.0/lib64
	TARGET_DEVICE := $(TARGET_GPU)
	GPU := GTX980
        NGPU := 1
endif



# thread number restrictions:
# BDy <= BDx
# BDx % warpSize == 0
# BDY * BDx > MAXWEI
# BDy <= warpSize           (always true)
# BDY * BDx <= 1024   (always true)

ifeq ($(GPU), GTX780)
	DEVFLAGS := -gencode arch=compute_35,code=sm_35 -maxrregcount=54 -Xptxas -dlcm=ca
#	DEVFLAGS := -Xptxas -dlcm=ca
	GD := 24
	BDy := 8
	BDx := 128
else ifeq ($(GPU), K20X)
	DEVFLAGS := -gencode arch=compute_35,code=sm_35 -maxrregcount=54 -Xptxas -dlcm=ca
#	DEVFLAGS := -Xptxas -dlcm=ca
	GD := 28
	BDy := 8
	BDx := 128
else ifeq ($(GPU), GTX980)
	DEVFLAGS := -gencode arch=compute_35,code=sm_52 -maxrregcount=61 -Xptxas -dlcm=ca
	GD := 32
	BDy := 8
	BDx := 64
else
	GD := 16
	BDy := 8
	BDx := 64
endif


# CUDA warp parameters 
WARP_SZ := 32

# Threads Per Block
TperB := $(shell echo $(BDx)*$(BDy) | bc)
WARPperB := $(shell echo $(TperB)/$(WARP_SZ) | bc)
WARPXperB := $(shell echo $(BDx)/$(WARP_SZ) | bc)
# rounding up to nearest power of 2
#BDx_POWER2 := $(shell echo "from math import log; a=2**int(log($(BDx),2)); print a" | python)
#TperB_POWER2 := $(shell echo "from math import log; a=2**int(log($(TperB),2)); print a" | python)

DMARCRO_GPU += -DGD=$(GD)
DMARCRO_GPU += -DBDx=$(BDx)
DMARCRO_GPU += -DBDy=$(BDy)
#DMARCRO_GPU += -DTperB=$(TperB)
DMARCRO_GPU += -DWARPperB=$(WARPperB)
DMARCRO_GPU += -DWARPXperB=$(WARPXperB)
#DMARCRO_GPU += -DBDx_POWER2=$(BDx_POWER2)
#DMARCRO_GPU += -DTperB_POWER2=$(TperB_POWER2)
DMARCRO_GPU += -DNGPU=$(NGPU)





EXE := dock
OBJ_HOST := main.o dock_cs.o load.o data.o rmsd.o util.o util_optimize_soa.o util_print.o
#OBJ_HOST := dock.o load.o data.o rmsd.o util.o util_optimize_aos.o util_print.o

DMARCRO := -DTARGET_GPU=$(TARGET_GPU) -DTARGET_CPU=$(TARGET_CPU) -DTARGET_MIC=$(TARGET_MIC) -DTARGET_DEVICE=$(TARGET_DEVICE)

HEADPATH += -I./module
HEADPATH += -I/usr/local/packages/cuda/6.5/include
LIBPATH += -L./module/yeah
LINKFLAGS += -lm -lyeahcpp

DEVFLAGS += $(DMARCRO) $(HEADPATH)
HOSTFLAGS += $(DMARCRO) $(HEADPATH)


MY_INTEL_BASE := /usr/local/compilers/Intel/cluster_studio_xe_2013.1.046/impi/4.1.3.048/intel64
MPIFLAGS += -I${MY_INTEL_BASE}/include
MPIFLAGS += -L${MY_INTEL_BASE}/lib
MPIFLAGS += -Xlinker --enable-new-dtags
MPIFLAGS += -Xlinker -rpath -Xlinker {MY_INTEL_BASE}/lib
MPIFLAGS += -lmpigc4 -lmpigf -lmpi -lmpigi -ldl -lrt -lpthread



ifeq ($(TARGET_DEVICE), $(TARGET_GPU))
#	CPP_HOST := icpc $(MPIFLAGS)
#	FC_HOST := ifort $(MPIFLAGS)

	CPP_HOST := g++
	FC_HOST := gfortran
	CPP_DEV := nvcc
	SRC_DEV := kernel_cuda.cu
	OBJ_DEV := kernel_cuda.o

	DEVFLAGS += -O3
	DEVFLAGS += $(DMARCRO_GPU)
	DEVFLAGS += -Xptxas -v
#	DEVFLAGS += #-Xcompiler -fopenmp
#	DEVFLAGS += -keep --keep-dir tmp -G --source-in-ptx
#
	HOSTFLAGS += -Wall -O3

	LINKFLAGS += -lcudart #-dlink

else ifeq ($(TARGET_DEVICE), $(TARGET_CPU))
#	CPP_HOST := g++
#	FC_HOST := gfortran
#	CPP_DEV := g++
#	SRC_DEV := kernel_cpp.C
#	OBJ_DEV := kernel_cpp.o
#	DEVFLAGS += -fopenmp
#	LINKFLAGS += -fopenmp


#	CPP_HOST := mpic++
#	FC_HOST := mpif90
#	CPP_DEV := mpic++


	CPP_HOST := icpc $(MPIFLAGS)
	FC_HOST := ifort $(MPIFLAGS)
	CPP_DEV := icpc $(MPIFLAGS)
	SRC_DEV := kernel_cpp.C
	OBJ_DEV := kernel_cpp.o

	DEVFLAGS += -Wall
	DEVFLAGS += -O3 #-ip -ipo -fma
#	DEVFLAGS += -xAVX
	DEVFLAGS += -fno-fnalias -ansi-alias -fargument-noalias
	DEVFLAGS += -openmp
	DEVFLAGS += -vec-threshold0
#	DEVFLAGS += -no-vec

	HOSTFLAGS += -Wall -O3

#        LINKFLAGS += -static
#	LINKFLAGS += -qopenmp
#	LINKFLAGS += -qopt-report=3
#	LINKFLAGS += -qopt-report-phase=vec
#	LINKFLAGS += -qopt-report-file=tmp/optrpt.txt

        LINKFLAGS += -static
	LINKFLAGS += -openmp
	LINKFLAGS += -opt-report=3
	LINKFLAGS += -opt-report-phase=vec
	LINKFLAGS += -opt-report-file=tmp/optrpt.txt



else ifeq ($(TARGET_DEVICE), $(TARGET_MIC))
	CPP_HOST := icc -mmic
	FC_HOST := ifort -mmic
	CPP_DEV := icc -mmic
	SRC_DEV := kernel_cpp.C
	OBJ_DEV := kernel_cpp.o
	LINKFLAGS += -openmp
	DEVFLAGS += -openmp
endif





default: $(EXE)

buildlib:
	cd modules/yeah && make clean
	cd modules/yeah && make

$(EXE): $(OBJ_HOST) $(OBJ_DEV)
	$(CPP_HOST) $(HOSTFLAGS) $(LIBPATH) $(OBJ_HOST) $(OBJ_DEV) -o $@ $(LINKFLAGS)

$(OBJ_DEV): $(SRC_DEV) kernel_cuda_l1_initcurand.cu kernel_cuda_l1_montecarlo.cu kernel_cuda_l2_util.cu kernel_cuda_prologue_epilogue.cu
	$(CPP_DEV) $(DEVFLAGS) -c $<

%.o: %.C
	$(CPP_HOST) $(HOSTFLAGS) -c $<

%.o: %.f
	$(FC_HOST) -O3 -c $<

#dock.o: dock.cu
#	$(CPP_DEV) $(DEVFLAGS) -c $<


install: $(EXE)
	(cp $(EXE) ../data/)
	(cp $(EXE) ../bin/)

run: $(EXE)
	cd ../data && ../src/$(EXE)

echo:
	@(echo BDx $(BDx))
	@(echo BDy $(BDy))
	@(echo TperB $(TperB))
	@(echo WARPperB $(WARPperB))
	@(echo WARPXperB $(WARPXperB))

clean:
	@(rm -f ${EXE} ${OBJ_HOST} ${OBJ_DEV} count *i *cubin* *.ptx *cudafe* *.fatbin* *.hash)






