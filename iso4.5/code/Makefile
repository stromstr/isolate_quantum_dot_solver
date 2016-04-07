LIBDIR  = -L/MKLPATH -L/opt/intel/lib/intel64/ -L/usr/lib64/ -L/usr/lib/ -I/MKLINCLUDE
LIBS    = $(LIBDIR)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
OBJS = spec.o calcams.o buildoperator.o swapspin.o swapms.o
FILE = .
NAME = ${HOME}/work/pro/iso4.5/hsms
NAME: ${OBJS}
	ifort  -o ${NAME} ${OBJS} ${LIBS}
spec.o: ${FILE}/spec.f90
	ifort  -c  ${FILE}/spec.f90
calcams.o: ${FILE}/calcams.f90
	ifort  -c  ${FILE}/calcams.f90
buildoperator.o: ${FILE}/buildoperator.f90
	ifort  -c ${FILE}/buildoperator.f90
swapspin.o: ${FILE}/swapspin.f90
	ifort  -c ${FILE}/swapspin.f90
swapms.o: ${FILE}/swapms.f90
	ifort  -c ${FILE}/swapms.f90
clean:
	rm -f ${NAME} ${OBJS}
