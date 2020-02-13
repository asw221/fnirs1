SRC = mcmc.cpp cholesky.cpp randgen.cpp mybspline.cpp hrf.cpp knots.cpp statistics.cpp config_info.cpp dlm.cpp kernel_reg.cpp
OBJ = mcmc.o cholesky.o randgen.o mybspline.o hrf.o knots.o statistics.o config_info.o dlm.o kernel_reg.o
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11
## INC = -I/Applications/MATLAB_R2019a.app/extern/include
INC = 
LFLAGS = 
LIBS = -lm -lfftw3
gp : ${OBJ} 
	$(CXX) -o foo ${OBJ} ${LIBS} ${INC}

${OBJ} : ${SRC}
	$(CXX) -c ${CXXFLAGS}  ${SRC} ${INC}

{OBJ} : cholesky.h randgen.h

clean:
	rm ${OBJ}
	rm *~
