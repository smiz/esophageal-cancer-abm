# Compiler flags
ADEVS=${HOME}/Code/adevs-code
ABM=${PWD}
CFLAGS = -O3 -fopenmp -Wall -std=c++11 -I${ADEVS}/include -I${ABM}
LIBS = \
	-lgsl \
	-lblas 

# Best bet for GNU compiler
CXX = g++
OBJS = \
	   common.o \
	   TissueVolume.o \
		main.o

.SUFFIXES: .cpp 
.cpp.o:
	${CXX} ${CFLAGS} ${OPTFLAG} -o $@ -c $<

objs: ${OBJS}
	${CXX} ${CFLAGS} ${OBJS} ${LIBS}

test: common.o
	${CXX} ${CFLAGS} test_common.cpp common.o ${LIBS}
	./a.out
	
clean:
	rm -f *.o a.out *csv*
