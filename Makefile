#############################################
### MAKE file for OMEGA3 Metagenome Assembler
#############################################

-include ../makefile.init


BUILD_GRAPH_CODE_DIR = src/BuildGraph/Release
BUILD_GRAPH_MPI_CODE_DIR = src/BuildGraphMPI/Release
BUILD_GRAPH_MPIRMA_CODE_DIR = src/BuildGraphMPIRMA/Release
SIMPLIFY_GRAPH_CODE_DIR = src/SimplifyGraph/Release

.PHONY: all clean 

all:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR)
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR)
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR)
		cp $(BUILD_GRAPH_CODE_DIR)/buildG .
		cp $(BUILD_GRAPH_MPI_CODE_DIR)/buildG-MPI .
		cp $(BUILD_GRAPH_MPI_CODE_DIR)/buildG-MPIRMA .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .
		mkdir Omega3
		cp buildG Omega3
		cp buildG-MPI Omega3
		cp buildG-MPIRMA Omega3
		cp fullsimplify Omega3
		cp parsimplify Omega3
		cp omega3.cfg Omega3
		cp runOmega3.sh Omega3
		cp README.md Omega3
		cp OUTPUT.md Omega3
		cp LICENSE Omega3
		tar -cvzf Omega3_x86-Linux.tar.gz Omega3
		rm -rf Omega3
clean:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR) clean
		$(MAKE) -C $(BUILD_GRAPH_MPI_CODE_DIR) clean
		$(MAKE) -C $(BUILD_GRAPH_MPIRMA_CODE_DIR) clean
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) clean
		-$(RM) -rf buildG buildG-MPI buildG-MPIRMA fullsimplify parsimplify Omega3_x86-Linux.tar.gz Omega3



