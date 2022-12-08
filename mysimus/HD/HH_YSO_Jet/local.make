LIBS += hdf5_serial
LIB_DIRS += $(LOCAL_LIBS)/hdf5/serial/lib 
INC_DIRS += $(LOCAL_LIBS)/hdf5/serial/include

LIBS += grackle
LIB_DIRS += /obs/mrabenanahary/MPI-AMRVAC/amrvac/src/local/lib
INC_DIRS += $(LD_INC_PATH)

LIBS77 += grackle
LIB_DIRS77 += /obs/mrabenanahary/MPI-AMRVAC/amrvac/src/local/lib
INC_DIRS77 += $(LD_INC_PATH)
