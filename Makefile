#############################################
### MAKE file for DISCO Metagenome Assembler
#############################################

-include ../makefile.init

ifeq ($(READGZ), 1)  #at this point, the makefile checks if FEATURE is enabled
GZOPTS = 1  #variable passed to g++
else
GZOPTS = 0  #variable passed to g++
endif


BUILD_GRAPH_CODE_DIR = src/BuildGraph/Release
BUILD_GRAPH_MPI_CODE_DIR = src/BuildGraphMPI/Release
BUILD_GRAPH_MPIRMA_CODE_DIR = src/BuildGraphMPIRMA/Release
SIMPLIFY_GRAPH_CODE_DIR = src/SimplifyGraph/Release

.PHONY: all clean 

openmp:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		cp $(BUILD_GRAPH_CODE_DIR)/buildG .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
mpi-dist-comp:
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
		cp $(BUILD_GRAPH_MPI_CODE_DIR)/buildG-MPI .
mpi-dist-mem:
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
		cp $(BUILD_GRAPH_MPIRMA_CODE_DIR)/buildG-MPIRMA .
all:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		cp $(BUILD_GRAPH_CODE_DIR)/buildG .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
		cp $(BUILD_GRAPH_MPI_CODE_DIR)/buildG-MPI .
		cp $(BUILD_GRAPH_MPIRMA_CODE_DIR)/buildG-MPIRMA .
release:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		cp $(BUILD_GRAPH_CODE_DIR)/buildG .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
		cp $(BUILD_GRAPH_MPI_CODE_DIR)/buildG-MPI .
		cp $(BUILD_GRAPH_MPIRMA_CODE_DIR)/buildG-MPIRMA .
		mkdir Disco
		cp buildG Disco
		cp buildG-MPI Disco
		cp buildG-MPIRMA Disco
		cp fullsimplify Disco
		cp parsimplify Disco
		cp disco.cfg Disco
		cp disco_2.cfg Disco
		cp disco_3.cfg Disco
		cp runDisco.sh Disco
		cp runAssembly.sh Disco
		cp runAssembly-MPI.sh Disco
		cp runECC.sh Disco
		cp runDisco-MPI.sh Disco
		cp runDisco-MPI-ALPS.sh Disco
		cp runDisco-MPI-SLURM.sh Disco
		cp README.md Disco
		cp OUTPUT.md Disco
		cp LICENSE Disco
		cp -r bbmap Disco
		tar -cvzf Disco_All_x86-Linux.tar.gz Disco
		rm -rf Disco
clean:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR) clean
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR) clean
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR) clean
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) clean
		-$(RM) -rf buildG buildG-MPI buildG-MPIRMA fullsimplify parsimplify 
		-$(RM) -rf Disco_omp_x86-Linux.tar.gz Disco_MPI_DC_x86-Linux.tar.gz 
		-$(RM) -rf Disco_All_x86-Linux.tar.gz Disco_MPI_DM_x86-Linux.tar.gz Disco