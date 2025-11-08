#######################################################################
# This is a simple Makefile for compiling your own channelflow programs.
# To compile and link program foo.cpp, run "make foo". It is best to do this
# outside the channelflow source tree. 
#
# You might need to modify the next few variables to specify the locations
# of header files and libraries on your system. In particular you will 
# want to change CHFLOWINC and CHFLOWLIB to point to your installation 
# of channelflow instead of mine (in ~gibson). 

# This configuration is for trillian, UNH's Cray XE6m-200, using HDF5 & Eigen
# installed in ~gibson/packages and FFTW installed in system directories
#CHFLOWINC = /home/math/gibson/channelflow/include
#CHFLOWLIB = /home/math/gibson/channelflow/lib
#EIGENINC = /home/math/gibson/packages/eigen-3.2.4/
#HDF5INC  = /home/math/gibson/packages/hdf5-1.8.15/include
#HDF5LIB  = /home/math/gibson/packages/hdf5-1.8.15/lib
#FFTWINC  = /opt/fftw/3.3.0.2/x86_64/include
#FFTWLIB  = /opt/fftw/3.3.0.2/x86_64/lib
#INCLUDES = -I$(CHFLOWINC) -I$(FFTWINC) -I$(EIGENINC) -I$(HDF5INC) 
#LIBPATHS = -L$(CHFLOWLIB) -L$(FFTWLIB) -L$(HDF5LIB) -Wl,-rpath,$(CHFLOWLIB):$(FFTWLIB):$(HDF5LIB)
LIBS = -lchflow -lfftw3 -lhdf5 -lhdf5_cpp -lz -lszip -lm

# This configuration is for typical Linux distributions. Here packages are installed 
# in predictable system locations that the compiler looks in automatically. 
EIGENINC  = /usr/include/eigen3
CHFLOWINC = /home/gibson/channelflow-2.0c/include
CHFLOWLIB = /home/gibson/channelflow-2.0c/lib
INCLUDES = -I$(CHFLOWINC) -I$(EIGENINC)
LIBPATHS = -L$(CHFLOWLIB) -Wl,-rpath,$(CHFLOWLIB)
LIBS = -lchflow -lfftw3 -lhdf5 -lhdf5_cpp -lz -lm



#######################################################################
# You should probably not have to modify anything below this line

# Define default C flags for optimized and debugging compilation
# Possible additional flags: -mtune=pentium4, -mtune-=athlon64
CFLAGS = -w -Wall -O3 -fomit-frame-pointer -DNDEBUG 

# Make optimized executables from objects and libraries
% : %.cpp 
	mpicxx $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBPATHS) $(LIBS)

# Make optimized objects from source
%.o : %.cpp 
	mpicxx $(CFLAGS) $(INCLUDES) -c $<

clean :
	rm -f *.o 
