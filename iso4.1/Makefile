LIBDIR  = -L/MKLPATH -L/opt/intel/lib/intel64/ -L/usr/lib64/ -L/usr/lib/ -I/MKLINCLUDE
LIBS    = $(LIBDIR)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
DIR_SRC = ./src
DIR_MOD = ./module
DIR_OBJECTS = ./objects
DIR_OBJMOD  = ./objmod
DIR_BIN = ./bin

SRC     = $(wildcard $(DIR_SRC)/*.f90)
SRC_MOD  = $(wildcard $(DIR_MOD)/*.f90)
OBJECTS = $(patsubst %.f90,$(DIR_OBJECTS)/%.o,$(notdir $(SRC)))
OBJMOD  = $(patsubst %.f90,$(DIR_OBJMOD)/%.mod,$(notdir $(SRC_MOD)))
TARGET  = showy
BIN_TARGET = ${DIR_BIN}/${TARGET}

F77     = ifort

${BIN_TARGET} : ${OBJMOD} ${OBJECTS}
	${F77} -module objmod objects/*.o ${LIBS} -o ${BIN_TARGET}
	echo compling done!

${DIR_OBJMOD}/%.mod: ${DIR_MOD}/%.f90
	${F77} -c -module objmod $<

${DIR_OBJECTS}/%.o : ${DIR_SRC}/%.f90
	${F77} -c -module objmod  $<
	mv ./*.o ${DIR_OBJECTS}

#${BIN_TARGET} : ${OBJECTS}
#	${F77} ${LIBS} ${OBJECTS} -o $@
#
clean : 
	rm -f ${BIN_TARGET} ${OBJECTS} ${OBJMOD}
