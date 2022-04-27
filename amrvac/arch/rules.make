VACPP := $(AMRVAC_DIR)/src/vacpp.pl

# Disable built-in make rules
.SUFFIXES: 

# How to compile modules into object files
mod_%.o: mod_%.f
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

# How to get .mod files from modules. Modules are automatically updated, and
# only need to be explicitly generated when they were manually removed.
mod_%.mod: mod_%.f mod_%.o
	@test -f $@ || $(F90) $(F90FLAGS) -c $(@:.mod=.f) -o $(@:.mod=.o) $(addprefix -I,$(INC_DIRS))

# How to generate object files
%.o: %.f
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

# How to translate .t source files to normal Fortran
%.f: %.t
	$(VACPP) $(PPFLAGS) -d=$(NDIM) $< > $(@)

# How to generate executables
%: %.o
	$(LINK) $(F90FLAGS) $^ -o $@ $(addprefix -L,$(LIB_DIRS)) $(addprefix -l,$(LIBS)) $(addprefix -L,$(LIB_DIRS77)) $(addprefix -l,$(LIBS77))
	
	
# How to compile modules into object files
mod_%chemistry.o: mod_%chemistry.f
	$(F77) $(F77FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS)) $(addprefix -I,$(INC_DIRS77))

# How to get .mod files from modules. Modules are automatically updated, and
# only need to be explicitly generated when they were manually removed.
mod_%chemistry.mod: mod_%chemistry.f mod_%chemistry.o
	@test -f $@ || $(F77) $(F77FLAGS) -c $(@:.mod=.f) -o $(@:.mod=.o) $(addprefix -I,$(INC_DIRS)) $(addprefix -I,$(INC_DIRS77))

# How to generate object files
%chemistry.o: %chemistry.f
	$(F77) $(F77FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS)) $(addprefix -I,$(INC_DIRS77))

# How to translate .t source files to normal Fortran
%chemistry.f: %chemistry.t
	$(VACPP) $(PPFLAGS) -d=$(NDIM) $< > $(@)

# How to generate executables
%: %chemistry.o
	$(LINK77) $(F77FLAGS) $^ -o $@ $(addprefix -L,$(LIB_DIRS)) $(addprefix -l,$(LIBS)) $(addprefix -L,$(LIB_DIRS77)) $(addprefix -l,$(LIBS77))	
	
	
	
# How to compile modules into object files
amrvac.o: amrvac.f
	$(F90_amrvac) $(F90FLAGS_amrvac) -c $< -o $@ $(addprefix -I,$(INC_DIRS)) $(addprefix -I,$(INC_DIRS77))

# How to get .mod files from modules. Modules are automatically updated, and
# only need to be explicitly generated when they were manually removed.
amrvac.mod: amrvac.f amrvac.o
	@test -f $@ || $(F90_amrvac) $(F90FLAGS_amrvac) -c $(@:.mod=.f) -o $(@:.mod=.o) $(addprefix -I,$(INC_DIRS)) $(addprefix -I,$(INC_DIRS77))

# How to translate .t source files to normal Fortran
amrvac.f: amrvac.t
	$(VACPP) $(PPFLAGS) -d=$(NDIM) $< > $(@)
	
# How to generate executables
amrvac: amrvac.o
	$(LINK_amrvac) $(F90FLAGS_amrvac) $^ -o $@ $(addprefix -L,$(LIB_DIRS)) $(addprefix -l,$(LIBS)) $(addprefix -L,$(LIB_DIRS77)) $(addprefix -l,$(LIBS77))
		
