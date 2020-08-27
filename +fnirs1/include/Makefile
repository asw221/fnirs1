
CXX      := clang++
SRC      := main.cpp mcmc.cpp cholesky.cpp randgen.cpp mybspline.cpp hrf.cpp kernel_reg.cpp knots.cpp statistics.cpp config_info.cpp dlm.cpp
OBJ      := main.o mcmc.o cholesky.o randgen.o mybspline.o hrf.o kernel_reg.o knots.o statistics.o config_info.o dlm.o
CXXFLAGS := -O3
INC      := -I/usr/local/include
LIB      := -L/usr/local/lib
LINK     := -lfftw3 -lfftw3_threads -lm
OMPLINK  := 

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LIB) $(LINK) -o fnirsdlm

$(OBJ): $(SRC)
	$(CXX) -c $(CXXFLAGS) $(SRC)

{OBJ}: cholesky.h randgen.h

clean:
	rm -f fnirsdlm *.o
