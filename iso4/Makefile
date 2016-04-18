LIBDIR  = -L/MKLPATH -L/opt/intel/lib/intel64/ -L/usr/lib64/ -L/usr/lib/ -I/MKLINCLUDE
LIBS    = $(LIBDIR)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_BIN = ./bin

SRC     = $(wildcard $(DIR_SRC)/*.f90)
OBJ     = $(patsubst %.f90,$(DIR_OBJ)/%.o,$(notdir $(SRC)))
TARGET  = hsms	
BIN_TARGET = ${DIR_BIN}/${TARGET}

F77     = ifort

${BIN_TARGET} : ${OBJ}
	${F77} ${LIBS} ${OBJ} -o $@

${DIR_OBJ}/%.o : ${DIR_SRC}/%.f90
	${F77} -c $<  -o $@

clean : 
	rm -f ${BIN_TARGET} ${OBJ} *.mod
