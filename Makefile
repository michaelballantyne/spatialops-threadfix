#
# USE THIS FOR LINUX:
#
TRILINOS_INCLUDE = /home/sutherland/apps/trilinos_jcs_opt_mpi/include
TRILINOS_LIB      = /home/sutherland/apps/trilinos_jcs_opt_mpi/lib 
LIBDIRS = -L./ -L$(TRILINOS_LIB)
EXTRA_LIBS = 

#
# USE THIS FOR MAC
#
#TRILINOS_INCLUDE = /jcs/software/trilinos/include
#TRILINOS_LIB     = /jcs/software/trilinos/lib
#LIBDIRS = -L./ -L$(TRILINOS_LIB) -L/sw/lib/gcc-lib/i386-apple-darwin8/4.0.3
#EXTRA_LIBS = -lf95

INCDIRS = -I./include -I$(TRILINOS_INCLUDE)

EPETRA_LIBS = -lepetra -lepetraext -lblas -llapack
AZTECOO_LIBS = -laztecoo -lteuchos 
LIBS = $(AZTECOO_LIBS)  $(EPETRA_LIBS) $(EXTRA_LIBS)

#CXXFLAGS = -O3 -Wall -fexpensive-optimizations -funroll-loops
CXXFLAGS = -g -Wall -DHAVE_MPI
#COMPILE_CXX = g++ -c $(CXXFLAGS) $(INCDIRS)
COMPILE_CXX = mpiCC -c $(CXXFLAGS) $(INCDIRS)

LINK = g++ $(CXXFLAGS) $(INCDIRS) $(LIBDIRS)

default: exe

OBJS =		\
	FVStaggeredSpatialOps.o \
	SpatialField.o \
	SpatialOperator.o \
	LinearSystem.o

FVStaggeredSpatialOps.o: ./src/FVStaggeredSpatialOps.cpp ./include/SpatialField.h ./include/SpatialOperator.h
	$(COMPILE_CXX) ./src/FVStaggeredSpatialOps.cpp

SpatialOperator.o: ./src/SpatialOperator.cpp ./include/SpatialOperator.h ./include/SpatialField.h
	$(COMPILE_CXX) ./src/SpatialOperator.cpp

SpatialField.o: ./src/SpatialField.cpp ./include/SpatialField.h ./include/SpatialOperator.h
	$(COMPILE_CXX) ./src/SpatialField.cpp

LinearSystem.o: ./src/LinearSystem.cpp ./include/SpatialOperator.h ./include/SpatialField.h ./include/LinearSystem.h
	$(COMPILE_CXX) ./src/LinearSystem.cpp




lib: $(OBJS)
	ar -r ./libspatialops.a $(OBJS)

exe: lib ./src/test.cpp
	$(LINK) ./src/test.cpp -lspatialops $(LIBS) -o test.x

clean: ; @rm *.o libspatialops.a test.x
