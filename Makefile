CC=mpicc
EBMS : EBMS.o runtime_parameters.o matrix.o
	${CC} -g -o EBMS EBMS.o runtime_parameters.o matrix.o
matrix.o : matrix.c matrix.h
	${CC} -g -c matrix.c -I.
runtime_parameters.o : runtime_parameters.c runtime_parameters.h
	${CC} -g -c runtime_parameters.c -I.
EBMS.o : EBMS.c EBMS.h
	${CC} -g -c EBMS.c -I.	
run : EBMS
	qsub -A Reactor_Modeling -q devel-ratio-16 -n 12 --proccount 12 -t 5 -O test1 ./EBMS params.in 
clean:
	\rm -f *.o *~ EBMS *output *error *cobaltlog