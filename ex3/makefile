OBJ    = A_times_u2.cc
CC     = icc
CFLAGS = -openmp -O0 -o
CLIBS  = -lmkl -lm -lgomp

A_mult_u.o: $(OBJ)
	$(CC) $(CFLAGS) A_times_u $(OBJ) $(CLIBS)

clean:
	rm -f *.o *~ A_times_u
