SRC = ./linalg.cpp ./ssh.cpp

HDR = $(SRC:.cpp=.hpp)
HDR += ./timing.hpp  

CC = g++ -std=c++14 -O2 -I./ 
CC += -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

CC += -fmax-errors=1 -fopenmp

LIB = -lm -larpack -lgfortran -lopenblas -lm
LIB +=  -lboost_program_options

CC += -I /opt/OpenBLAS/include/ -L /opt/OpenBLAS/lib/

spectral: spectral.cpp $(SRC) $(HDR)
	$(CC) $(SRC) ./$< $(LIB) -o ./$@
	
clean:
	rm -f -v ./spectral