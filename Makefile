#
# USE THIS FOR LINUX:
#
BOOST_INCLUDE = /home/sutherland/packages/boost_1_35_0
TRILINOS_INCLUDE = /home/sutherland/apps/trilinos_jcs_opt/include
TRILINOS_LIB      = /home/sutherland/apps/trilinos_jcs_opt/lib 
DAIXT_INCLUDE     = /home/sutherland/apps/daixtrose/include
LIBDIRS = -L./ -L$(TRILINOS_LIB)
EXTRA_LIBS = 
CXX = g++4

#
# USE THIS FOR MAC
#
#BOOST_INCLUDE = /jcs/software/boost_1_35_0
#TRILINOS_INCLUDE = /jcs/software/trilinos/include
#TRILINOS_LIB     = /jcs/software/trilinos/lib
#DAIXT_INCLUDE    = /jcs/software/daixtrose-0.0.3/jcs_install/include
#LIBDIRS = -L./ -L$(TRILINOS_LIB) -L/sw/lib/gcc-lib/i386-apple-darwin8/4.0.3
#EXTRA_LIBS = -lf95
#CXX = g++

INCDIRS = -I./include -I$(DAIXT_INCLUDE) -I$(TRILINOS_INCLUDE) -I$(BOOST_INCLUDE)

EPETRA_LIBS = -lepetra -lepetraext -lblas -llapack
AZTECOO_LIBS = -laztecoo -lteuchos 
LIBS = $(AZTECOO_LIBS)  $(EPETRA_LIBS) $(EXTRA_LIBS)

CXXFLAGS = -O4 -Wall -fexpensive-optimizations -funroll-loops -DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR -DNDEBUG -DUINTAH_FIELD_TYPES #-DSAMRAI_FIELD_TYPES
#CXXFLAGS = -g -Wall -O0 -DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR -DUINTAH_FIELD_TYPES #-DSAMRAI_FIELD_TYPES 
COMPILE_CXX = $(CXX) -c $(CXXFLAGS) $(INCDIRS)
#COMPILE_CXX = mpiCC -c $(CXXFLAGS) $(INCDIRS)

LINK = $(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS)

default: testnew poisson simplebc
all: testnew poisson simplebc

OBJS =		\
	LinAlgTrilinos.o \
	LinAlgUBlas.o \
	LinearSystem.o

bcTest.o: ./test/bcTest.cpp ./include/*.h
	$(COMPILE_CXX) -I ./ ./test/bcTest.cpp

buildOps.o: ./test/buildOps.cpp ./include/*.h
	$(COMPILE_CXX) -I ./test ./test/buildOps.cpp

testNew.o: ./test/testNew.cpp ./include/*.h ./test/*.h
	$(COMPILE_CXX) -I./test ./test/testNew.cpp

LinAlgTrilinos.o: ./src/LinAlgTrilinos.cpp ./include/LinAlgTrilinos.h
	$(COMPILE_CXX) ./src/LinAlgTrilinos.cpp

LinAlgUBlas.o: ./src/LinAlgUBlas.cpp ./include/LinAlgUBlas.h
	$(COMPILE_CXX) ./src/LinAlgUBlas.cpp

LinearSystem.o: ./src/LinearSystem.cpp ./include/SpatialOperator.h ./include/SpatialField.h ./include/LinearSystem.h
	$(COMPILE_CXX) ./src/LinearSystem.cpp

testPoisson.o: ./test/testPoisson.cpp ./include/*.h ./test/*.h
	$(COMPILE_CXX) -I./test ./test/testPoisson.cpp

lib: $(OBJS)
	ar -r ./libspatialops.a $(OBJS); ranlib ./libspatialops.a

simplebc: lib bcTest.o
	$(LINK) bcTest.o -lspatialops $(LIBS) -o simplebc.x

testnew: lib testNew.o buildOps.o
	$(LINK) testNew.o buildOps.o -lspatialops $(LIBS) -o testnew.x

poisson: lib testPoisson.o buildOps.o
	$(LINK) testPoisson.o buildOps.o -lspatialops $(LIBS) -o testpoisson.x

clean: ; @rm *.o libspatialops.a testnew.x testpoisson.x
