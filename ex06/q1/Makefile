#CC = gcc
CC = mpiCC
CFLAGS = -I/cluster/work/infk/psumit/trilinos-9.0.1/LINUX_MPI/include
LDFLAGS = \
	-L/cluster/work/infk/psumit/ParMetis-3.1.1 \
        -L/cluster/work/infk/psumit/trilinos-9.0.1/LINUX_MPI/lib -lmetis  -lepetraext -ltriutils -lepetra -lteuchos  \
        -L ${MKL_HOME}/lib/em64t -lmkl_blacs_openmpi_lp64 -lmkl_em64t -lguide -lml

mat_vec_trilinos: mat_vec_trilinos.o 
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)


mat_vec_trilinos.o: mat_vec_trilinos.cpp
	$(CC) $(CFLAGS) -c mat_vec_trilinos.cpp



clean:
	rm *.o mat_vec_trilinos mat_vec_trilinos.o
