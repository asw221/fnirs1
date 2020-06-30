
CXX      := clang++
SRC      := main.cpp mcmc.cpp cholesky.cpp randgen.cpp mybspline.cpp hrf.cpp knots.cpp statistics.cpp config_info.cpp dlm.cpp kernel_reg.cpp marg_eta_updates.cpp
OBJ      := main.o mcmc.o cholesky.o randgen.o mybspline.o hrf.o knots.o statistics.o config_info.o dlm.o kernel_reg.o marg_eta_updates.o
CXXFLAGS := -g -O2
INC      := -I/usr/local/include
LIB      := -L/usr/local/lib
LINK     := -lfftw3 -lfftw3_threads -lm
OMPLINK  := -Xpreprocessor -fopenmp -lomp

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LIB) $(LINK) $(OMPLINK) -o fnirsdlm

$(OBJ): $(SRC)
	$(CXX) -c $(CXXFLAGS) $(SRC) $(OMPLINK)

{OBJ}: cholesky.h randgen.h

clean:
	rm -f fnirsdlm *.o
