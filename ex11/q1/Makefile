CPP = mpicxx
CFLAGS = 
LDFLAGS = 

q1: q1.o 
		$(CPP) -o $@ $@.o $(OBJECTS) $(LDFLAGS)


q1.o: q1.cpp
		$(CPP) $(CFLAGS) -c q1.cpp



clean:
		rm *.o q1 q1.o

