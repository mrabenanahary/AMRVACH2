ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

#openmpi
ARCH =default
NDIM := 2

# By exporting these can be used when building libamrvac
export ARCH NDIM
OBJ_DIR := Objects
SRC_DIR := $(AMRVAC_DIR)/src
LIB_DIR := $(AMRVAC_DIR)/lib/$(NDIM)d_$(ARCH)
LIB_MAKE := $(AMRVAC_DIR)/arch/lib.make
LIB_AMRVAC := $(LIB_DIR)/libamrvac.a



SRC_DIR+= $(OBJ_DIR)
# These are used for compilation
INC_DIRS := $(LIB_DIR)
LIB_DIRS := $(LIB_DIR)
LIBS := amrvac

.PHONY: all clean allclean force

all: amrvac

# Include architecture and compilation rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# Optionally include a local user makefile
-include local.make

# Where to find amrvac.t
vpath %.t $(SRC_DIR)

# Keep mod_usr.f for inspection
.PRECIOUS: mod_usr_yso_jet.f \
	        $(OBJ_DIR)/mod_obj_global_parameters.f $(OBJ_DIR)/mod_obj_mat.f \
					$(OBJ_DIR)/mod_obj_relativistic_wind.f $(OBJ_DIR)/mod_obj_star_envelope.f\
					$(OBJ_DIR)/mod_obj_star.f $(OBJ_DIR)/mod_obj_pulsar.f\
					$(OBJ_DIR)/mod_obj_sn_remnant.f $(OBJ_DIR)/mod_obj_ism.f \
					$(OBJ_DIR)/mod_obj_dust.f $(OBJ_DIR)/mod_obj_cloud.f  \
	 $(OBJ_DIR)/mod_obj_usr_unit.f

# Intermediate files are removed
.INTERMEDIATE: amrvac.o mod_usr_yso_jet.o  mod_usr_yso_jet.mod

# Always try to build/update the amrvac library
$(LIB_AMRVAC): force
	@mkdir -p $(LIB_DIR)
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE)

clean:
	@echo $(SRC_DIR)
	@echo 'Cleaning local objects ("make allclean" cleans libamrvac)'
	$(RM) amrvac amrvac.o mod_usr_supernovae_dust.o \
	      $(OBJ_DIR)/mod_obj_global_parameters.o  $(OBJ_DIR)/mod_obj_mat.o\
              $(OBJ_DIR)/mod_obj_relativistic_wind.o $(OBJ_DIR)/mod_obj_star_envelope.o\
              $(OBJ_DIR)/mod_obj_sn_remnant.o $(OBJ_DIR)/mod_obj_ism.o \
	      $(OBJ_DIR)/mod_obj_dust.o $(OBJ_DIR)/mod_obj_cloud.o\
	      $(OBJ_DIR)/mod_obj_star.o $(OBJ_DIR)/mod_obj_pulsar.o\
	      $(OBJ_DIR)/mod_obj_usr_unit.o\
	      $(OBJ_DIR)/mod_obj_global_parameters.mod  $(OBJ_DIR)/mod_obj_mat.mod\
				$(OBJ_DIR)/mod_obj_relativistic_wind.mod\
				$(OBJ_DIR)/mod_obj_star.mod $(OBJ_DIR)/mod_obj_pulsar.mod\
				$(OBJ_DIR)/mod_obj_sn_remnant.mod $(OBJ_DIR)/mod_obj_ism.mod \
				$(OBJ_DIR)/mod_obj_dust.mod $(OBJ_DIR)/mod_obj_cloud.mod\
				$(OBJ_DIR)/mod_obj_usr_unit.mod                         \
				mod_usr.mod \
				$(OBJ_DIR)/*.f mod_usr_yso_jet.f

# Also clean the library
allclean: clean
	@echo 'Cleaning libamrvac'
	@mkdir -p $(LIB_DIR)	# Prevent error message
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE) clean

# Dependencies
amrvac: $(OBJ_DIR)/mod_obj_global_parameters.o $(OBJ_DIR)/mod_obj_mat.o\
	      $(OBJ_DIR)/mod_obj_dust.o $(OBJ_DIR)/mod_obj_sn_remnant.o \
	      $(OBJ_DIR)/mod_obj_relativistic_wind.o $(OBJ_DIR)/mod_obj_star_envelope.o\
	      $(OBJ_DIR)/mod_obj_star.o $(OBJ_DIR)/mod_obj_pulsar.o \
	      $(OBJ_DIR)/mod_obj_ism.o $(OBJ_DIR)/mod_obj_cloud.o \
	      $(OBJ_DIR)/mod_obj_usr_unit.o \
              mod_usr_supernovae_dust.o amrvac.o

amrvac.o  mod_usr_yso_jet.o: $(LIB_AMRVAC)
$(OBJ_DIR)/mod_obj_dust.o: $(OBJ_DIR)/mod_obj_global_parameters.o $(OBJ_DIR)/mod_obj_mat.o \
	                   $(OBJ_DIR)/mod_obj_usr_unit.o
$(OBJ_DIR)/mod_obj_star_envelope.o:$(OBJ_DIR)/mod_obj_global_parameters.o $(OBJ_DIR)/mod_obj_mat.o \
	                           $(OBJ_DIR)/mod_obj_usr_unit.o
$(OBJ_DIR)/mod_obj_relativistic_wind.o: $(OBJ_DIR)/mod_obj_global_parameters.o $(OBJ_DIR)/mod_obj_mat.o \
	                                $(OBJ_DIR)/mod_obj_usr_unit.o
$(OBJ_DIR)/mod_obj_sn_remnant.o: $(OBJ_DIR)/mod_obj_global_parameters.o $(OBJ_DIR)/mod_obj_mat.o \
	                         $(OBJ_DIR)/mod_obj_usr_unit.o
$(OBJ_DIR)/mod_obj_star.o:$(OBJ_DIR)/mod_obj_global_parameters.o $(OBJ_DIR)/mod_obj_mat.o\
	                  $(OBJ_DIR)/mod_obj_usr_unit.o
$(OBJ_DIR)/mod_obj_pulsar.o:$(OBJ_DIR)/mod_obj_sn_remnant.o  $(OBJ_DIR)/mod_obj_relativistic_wind.o \
	                          $(OBJ_DIR)/mod_obj_star_envelope.o $(OBJ_DIR)/mod_obj_star.o\
         			  $(OBJ_DIR)/mod_obj_usr_unit.o
$(OBJ_DIR)/mod_obj_ism.o : $(OBJ_DIR)/mod_obj_dust.o
$(OBJ_DIR)/mod_obj_cloud.o:$(OBJ_DIR)/mod_obj_dust.o
mod_usr_yso_jet.o: $(OBJ_DIR)/mod_obj_global_parameters.o $(OBJ_DIR)/mod_obj_mat.o \
	$(OBJ_DIR)/mod_obj_dust.o $(OBJ_DIR)/mod_obj_relativistic_wind.o\
	$(OBJ_DIR)/mod_obj_sn_remnant.o $(OBJ_DIR)/mod_obj_ism.o $(OBJ_DIR)/mod_obj_cloud.o
amrvac.o: mod_usr_yso_jet.o
