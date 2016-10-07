#############################################
### MAKE file for OMEGA3 Metagenome Assembler
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
		mkdir Omega3
		cp buildG Omega3
		cp fullsimplify Omega3
		cp parsimplify Omega3
		cp omega3.cfg Omega3
		cp runOmega3.sh Omega3
		cp README.md Omega3
		cp OUTPUT.md Omega3
		cp LICENSE Omega3
		tar -cvzf Omega3_omp_x86-Linux.tar.gz Omega3
		rm -rf Omega3
mpi-dist-comp:
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
		cp $(BUILD_GRAPH_MPI_CODE_DIR)/buildG-MPI .
		mkdir Omega3
		cp buildG-MPI Omega3
		cp fullsimplify Omega3
		cp parsimplify Omega3
		cp omega3.cfg Omega3
		cp runOmega3-MPI.sh Omega3
		cp runOmega3-MPI-ALPS.sh Omega3
		cp runOmega3-MPI-SLRUM.sh Omega3
		cp README.md Omega3
		cp OUTPUT.md Omega3
		cp LICENSE Omega3
		tar -cvzf Omega3_MPI_DC_x86-Linux.tar.gz Omega3
		rm -rf Omega3
mpi-dist-mem:
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR) READGZ=$(GZOPTS)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) READGZ=$(GZOPTS)
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
		cp $(BUILD_GRAPH_MPIRMA_CODE_DIR)/buildG-MPIRMA .
		mkdir Omega3
		cp buildG-MPIRMA Omega3
		cp fullsimplify Omega3
		cp parsimplify Omega3
		cp omega3.cfg Omega3
		cp runOmega3-MPI.sh Omega3
		cp runOmega3-MPI-ALPS.sh Omega3
		cp runOmega3-MPI-SLRUM.sh Omega3
		cp README.md Omega3
		cp OUTPUT.md Omega3
		cp LICENSE Omega3
		tar -cvzf Omega3_MPI_DM_x86-Linux.tar.gz Omega3
		rm -rf Omega3
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
		mkdir Omega3
		cp buildG Omega3
		cp buildG-MPI Omega3
		cp buildG-MPIRMA Omega3
		cp fullsimplify Omega3
		cp parsimplify Omega3
		cp omega3.cfg Omega3
		cp runOmega3.sh Omega3
		cp runOmega3-MPI.sh Omega3
		cp runOmega3-MPI-ALPS.sh Omega3
		cp runOmega3-MPI-SLURM.sh Omega3
		cp README.md Omega3
		cp OUTPUT.md Omega3
		cp LICENSE Omega3
		tar -cvzf Omega3_All_x86-Linux.tar.gz Omega3
		rm -rf Omega3
clean:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR) clean
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR) clean
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR) clean
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) clean
		-$(RM) -rf buildG buildG-MPI buildG-MPIRMA fullsimplify parsimplify 
		-$(RM) -rf Omega3_omp_x86-Linux.tar.gz Omega3_MPI_DC_x86-Linux.tar.gz 
		-$(RM) -rf Omega3_All_x86-Linux.tar.gz Omega3_MPI_DM_x86-Linux.tar.gz Omega3