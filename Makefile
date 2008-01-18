#
# USE THIS FOR LINUX:
#
BOOST_INCLUDE = /home/sutherland/packages/boost_1_34_0
TRILINOS_INCLUDE = /home/sutherland/apps/trilinos_jcs_opt/include
TRILINOS_LIB      = /home/sutherland/apps/trilinos_jcs_opt/lib 
DAIXT_INCLUDE     = /home/sutherland/apps/daixtrose/include
LIBDIRS = -L./ -L$(TRILINOS_LIB)
EXTRA_LIBS = 

#
# USE THIS FOR MAC
#
#BOOST_INCLUDE = /jcs/software/boost_1_34_1
#TRILINOS_INCLUDE = /jcs/software/trilinos/include
#TRILINOS_LIB     = /jcs/software/trilinos/lib
#DAIXT_INCLUDE    = /jcs/software/daixtrose-0.0.3/jcs_install/include
#LIBDIRS = -L./ -L$(TRILINOS_LIB) -L/sw/lib/gcc-lib/i386-apple-darwin8/4.0.3
#EXTRA_LIBS = -lf95

INCDIRS = -I./include -I$(DAIXT_INCLUDE) -I$(TRILINOS_INCLUDE) -I$(BOOST_INCLUDE)

EPETRA_LIBS = -lepetra -lepetraext -lblas -llapack
AZTECOO_LIBS = -laztecoo -lteuchos 
LIBS = $(AZTECOO_LIBS)  $(EPETRA_LIBS) $(EXTRA_LIBS)

CXXFLAGS = -O4 -Wall -fexpensive-optimizations -funroll-loops -DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR -DNDEBUG -pg
#CXXFLAGS = -g -Wall -O0 -DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR
COMPILE_CXX = g++ -c $(CXXFLAGS) $(INCDIRS)
#COMPILE_CXX = mpiCC -c $(CXXFLAGS) $(INCDIRS)

LINK = g++ $(CXXFLAGS) $(INCDIRS) $(LIBDIRS)

default: testnew poisson
all: testnew poisson

OBJS =		\
	LinAlgTrilinos.o \
	LinAlgUBlas.o \
	LinearSystem.o

buildOps.o: ./src/test/buildOps.cpp ./include/*.h
	$(COMPILE_CXX) -I ./src/test ./src/test/buildOps.cpp

testNew.o: ./src/test/testNew.cpp ./include/*.h ./src/test/*.h
	$(COMPILE_CXX) -I./src/test ./src/test/testNew.cpp

LinAlgTrilinos.o: ./src/LinAlgTrilinos.cpp ./include/LinAlgTrilinos.h
	$(COMPILE_CXX) ./src/LinAlgTrilinos.cpp

LinAlgUBlas.o: ./src/LinAlgUBlas.cpp ./include/LinAlgUBlas.h
	$(COMPILE_CXX) ./src/LinAlgUBlas.cpp

LinearSystem.o: ./src/LinearSystem.cpp ./include/SpatialOperator.h ./include/SpatialField.h ./include/LinearSystem.h
	$(COMPILE_CXX) ./src/LinearSystem.cpp

test_alberto.o: ./src/test/test_alberto.cpp 
	$(COMPILE_CXX) ./src/test/test_alberto.cpp -I ./src/test/test_alberto

testPoisson.o: ./src/test/testPoisson.cpp ./include/*.h ./src/test/*.h
	$(COMPILE_CXX) -I./src/test ./src/test/testPoisson.cpp

lib: $(OBJS)
	ar -r ./libspatialops.a $(OBJS)

testnew: lib testNew.o buildOps.o
	$(LINK) testNew.o buildOps.o -lspatialops $(LIBS) -o testnew.x

poisson: lib testPoisson.o buildOps.o
	$(LINK) testPoisson.o buildOps.o -lspatialops $(LIBS) -o testpoisson.x

clean: ; @rm *.o libspatialops.a test.x
