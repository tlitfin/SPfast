include Makefile.in

SRC = protein2.cc align1.cc quatfit_theo.cc SPfast.cc quatfit.cc #SPfast
SRC0 = protein2.cc prepare_bin.cc #prepare_bin
SRC1 = protein2.cc extract_bin.cc #extract_bin
#SRC1 = protein2.cc align1.cc quatfit_theo.cc quatfit.cc

#
#pb11_includes := $(shell python3 -m pybind11 --includes)
#pb11_suffix := $(shell python3-config --extension-suffix)

all: SPfast.$(TYPE) prepare_bin.$(TYPE) extract_bin.$(TYPE) #SPlib$(pb11_suffix)
intel: all
gdb: all
gnu: all
mp: all
pg: all

# make object files with *.h dependencies
$(OBJ): *.h
$(OBJ0): *.h
$(OBJ1): *.h

CPPFLAGS = -fPIC -ffast-math -pthread -fopenmp -DMP
#CPPFLAGS1 = -std=gnu++11 -O3 -shared -fopenmp -DMP 

# link binaries to relevant object files
prepare_bin.$(TYPE): $(OBJ0)
	$(CPP) $(CPPFLAGS) $^ -o $@
extract_bin.$(TYPE): $(OBJ1)
	$(CPP) $(CPPFLAGS) $^ -o $@
SPfast.$(TYPE): $(OBJ)
	$(CPP) $(CPPFLAGS) $^ -o $@
#SPlib$(pb11_suffix): $(OBJ)
#	$(CPP) $(CPPFLAGS1) $^ -fPIC $(pb11_includes) SPlib.cpp -o $@
clean:
	rm -f .*.o
