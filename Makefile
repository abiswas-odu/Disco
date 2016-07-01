#############################################
### MAKE file for OMEGA3 Metagenome Assembler
#############################################

-include ../makefile.init


BUILD_GRAPH_CODE_DIR = src/BuildGraph/Release
SIMPLIFY_GRAPH_CODE_DIR = src/SimplifyGraph/Release

.PHONY: all clean 

all:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR)
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR)
		cp $(BUILD_GRAPH_CODE_DIR)/buildG .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/fullsimplify .
		cp $(SIMPLIFY_GRAPH_CODE_DIR)/parsimplify .

clean:
		$(MAKE) -C $(BUILD_GRAPH_CODE_DIR) clean
		$(MAKE) -C $(SIMPLIFY_GRAPH_CODE_DIR) clean
		-$(RM) buildG fullsimplify parsimplify



