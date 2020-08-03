
CXX      := clang++
SRC      := main.cpp mcmc.cpp cholesky.cpp randgen.cpp mybspline.cpp hrf.cpp knots.cpp statistics.cpp config_info.cpp dlm.cpp
OBJ      := main.o mcmc.o cholesky.o randgen.o mybspline.o hrf.o knots.o statistics.o config_info.o dlm.o
CXXFLAGS := -O3
INC      := -I/usr/local/include
LIB      := -L/usr/local/lib
LINK     := -lfftw3 -lfftw3_threads -lm -lomp
OMPLINK  := -Xpreprocessor -fopenmp

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LIB) $(LINK) $(OMPLINK) -o fnirsdlm

$(OBJ): $(SRC)
	$(CXX) -c $(CXXFLAGS) $(SRC) $(INC) $(OMPLINK)

{OBJ}: cholesky.h randgen.h

clean:
	rm -f fnirsdlm *.o
