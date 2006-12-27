INCDIRS = -I./include -I/home/sutherland/apps/trilinos/include
LIBDIRS = -L/jcs/software/trilinos/lib -L./ -L/home/sutherland/apps/trilinos/lib 
EPETRA_LIBS = -lepetra -lepetraext -lblas -llapack
AZTECOO_LIBS = -laztecoo -lteuchos 
LIBS =  $(EPETRA_LIBS) $(AZTECOO_LIBS)

CXXFLAGS = -g -Wall
COMPILE_CXX = g++ -c $(CXXFLAGS) $(INCDIRS)

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
	ar -r -o ./libspatialops.a $(OBJS)

exe: lib ./src/test.cpp
	$(LINK) ./src/test.cpp -lspatialops $(LIBS) -o test.x

clean: ; @rm *.o libspatialops.a test.x
