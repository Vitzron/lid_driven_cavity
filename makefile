CXX = g++
CXXFLAGS = -Wall -Ofast
LIBS = -I /usr/local/include/eigen3/Eigen		# change this path to Eigen directory
PROM = main
OBJS = U_BC.o V_BC.o P_BC.o Grid.o Grid_Utilities.o Eigen_Utilities.o\
	WENO.o Pressure.o main.o

${PROM}: ${OBJS}
	${CXX} -o ${PROM} ${OBJS} ${LIBS}
clean:
	rm -f *.o ${PROM}