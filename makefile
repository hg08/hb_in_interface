#=================================
#Do we want debug options for gcc?
#=================================
DEBUG='warning'
ifeq ($(DEBUG),'warning')
	CFLAGS= -Wall -Wextra -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wintrinsics-std -Wreal-q-constant -Wsurprising -Wunderflow -Wintrinsic-shadow -Walign-commons -ffrontend-optimize -pedantic -fcheck=all -g -pg
else ifeq ($(DEBUG),'fprof')
	CFLAGS= -Wall -Wextra -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wintrinsics-std -Wreal-q-constant -Wsurprising -Wunderflow -Wintrinsic-shadow -Walign-commons -ffrontend-optimize -pedantic -fcheck=all -g -pg -O0 -ftest-coverage -fprofile-arcs  #Compiler flags for debug run
else ifeq ($(DEBUG),'opt')
	CFLAGS= -O2
else
	CFLAGS=    #Compiler flags for normal run
endif

#=========
#Variables
#=========
SHELL=bash
FC=gfortran
SRC_DIR=src
OBJ_DIR=obj
EXEC=ihb_sulpizi
SRC_DUP= $(wildcard $(addprefix $(SRC_DIR)/, module.f95 *.f95)) #List of the sources (module.f95 has to be first) but it will be in duplicate
SRC=$(shell echo $(SRC_DUP) | tr ' ' '\n' | awk '!a[$$0]++' | tr '\n' ' ') #Deletion of the duplicate module.f95
OBJ= $(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.f95=.o))) #List of the object files 

#=================
#gfortran commands
#=================
# @$(command) here it is @$(FC) = silent mode. No line is written if everything is ok.
#Here ‘$<’ is the automatic variable that holds the name of the prerequisite; 
#‘$@’ is the automatic variable that holds the name of the target;
# See https://www.gnu.org/software/make/manual/make.html#Reading 

#Creation of the execuatble
$(EXEC): $(OBJ) #If target (main) is older than $(OBJ)
	@$(FC) -o $@ $(OBJ) $(CFLAGS) $(LDFLAGS)

#Creation of the objects
#Produce .o files if .f95 changed
#%.o: %.c
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f95
	@$(FC) -c -o $@ $< $(CFLAGS)


#============
#============


#Rebuild the dependances (even if another file is name clean)
.PHONY: clean 

#Remove the intermediate files
clean: 
	rm -rf $(OBJ_DIR)/*

#Remove the intermediate files
very_clean: 
	rm -rf $(OBJ_DIR)/* $(SRC_DIR)/TAGS mymod.mod vvaf

#Create tag files to navigate with emacs
#module.f95 is put at the end because usually, it is the one that I use the less
tag:
	etags $(SRC_DIR)/*.f95 $(SRC_DIR)/module.f95 -o $(SRC_DIR)/TAGS

#==================
#Debugging commands
#==================
#gprof
#~/bin/__sources_f90/__vvaf/vvaf input_vvaf
#gprof ~/bin/__sources_f90/__vvaf/vvaf | gprof2dot.py -n0 -e0 | dot -Tpng -o gprof_new.png
#gprof ~/bin/__sources_f90/__vvaf/vvaf > gprof_new.output

#gcov + lcov
#for i in *f95; do gcov $i ; done
#lcov --capture --directory ./ --output-file coverage.info
#genhtml coverage.info --output-directory out
