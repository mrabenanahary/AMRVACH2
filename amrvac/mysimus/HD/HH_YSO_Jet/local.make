LIBS += hdf5_serial
LIB_DIRS += $(LOCAL_LIBS)/hdf5/serial/lib 
INC_DIRS += $(LOCAL_LIBS)/hdf5/serial/include

LIBS += grackle
LIB_DIRS += $(LD_LIBRARY_PATH)
INC_DIRS += $(LD_INC_PATH)

LIBS77 += grackle
LIB_DIRS77 += $(LD_LIBRARY_PATH)
INC_DIRS77 += $(LD_INC_PATH)
