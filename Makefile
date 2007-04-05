#
# USE THIS FOR LINUX:
#
BOOST_INCLUDE = /home/sutherland/packages/boost_1_33_1
TRILINOS_INCLUDE = /home/sutherland/apps/trilinos_jcs_opt/include
TRILINOS_LIB      = /home/sutherland/apps/trilinos_jcs_opt/lib 
LIBDIRS = -L./ -L$(TRILINOS_LIB)
EXTRA_LIBS = 

#
# USE THIS FOR MAC
#
#BOOST_INCLUDE = /jcs/software/boost_1_33_1
#TRILINOS_INCLUDE = /jcs/software/trilinos/include
#TRILINOS_LIB     = /jcs/software/trilinos/lib
#LIBDIRS = -L./ -L$(TRILINOS_LIB) -L/sw/lib/gcc-lib/i386-apple-darwin8/4.0.3
#EXTRA_LIBS = -lf95

INCDIRS = -I./include -I$(TRILINOS_INCLUDE) -I$(BOOST_INCLUDE)

EPETRA_LIBS = -lepetra -lepetraext -lblas -llapack
AZTECOO_LIBS = -laztecoo -lteuchos 
LIBS = $(AZTECOO_LIBS)  $(EPETRA_LIBS) $(EXTRA_LIBS)

CXXFLAGS = -O3 -Wall -fexpensive-optimizations -funroll-loops
#CXXFLAGS = -g -Wall
COMPILE_CXX = g++ -c $(CXXFLAGS) $(INCDIRS)
#COMPILE_CXX = mpiCC -c $(CXXFLAGS) $(INCDIRS)

LINK = g++ $(CXXFLAGS) $(INCDIRS) $(LIBDIRS)

default: test

OBJS =		\
	FVStaggeredSpatialOps.o \
	LinAlgTrilinos.o \
	LinearSystem.o

test.o: ./src/test/test.cpp ./include/SpatialField.h ./include/SpatialOperator.h ./include/FVStaggeredSpatialOps.h ./include/LinAlgTrilinos.h ./include/FV2ndOrderTypes.h
	$(COMPILE_CXX) ./src/test/test.cpp

LinAlgTrilinos.o: ./src/LinAlgTrilinos.cpp ./include/LinAlgTrilinos.h
	$(COMPILE_CXX) ./src/LinAlgTrilinos.cpp

FVStaggeredSpatialOps.o: ./src/FVStaggeredSpatialOps.cpp ./include/SpatialField.h ./include/SpatialOperator.h
	$(COMPILE_CXX) ./src/FVStaggeredSpatialOps.cpp

LinearSystem.o: ./src/LinearSystem.cpp ./include/SpatialOperator.h ./include/SpatialField.h ./include/LinearSystem.h
	$(COMPILE_CXX) ./src/LinearSystem.cpp

test_alberto.o: ./src/test/test_alberto.cpp 
	$(COMPILE_CXX) ./src/test/test_alberto.cpp -I ./src/test/test_alberto


lib: $(OBJS)
	ar -r ./libspatialops.a $(OBJS)

test: lib test_alberto.o
	$(LINK) test_alberto.o  -lspatialops $(LIBS) -o test.x

smalltest: lib test.o
	$(LINK) test.o -lspatialops $(LIBS) -o smalltest.x

clean: ; @rm *.o libspatialops.a test.x
