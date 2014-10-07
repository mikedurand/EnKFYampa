FC=mpif90
FLAGS=-g -C

filter: *.f *.f90
	$(FC) $(FLAGS) -o $@ $^

clean:
	rm *.o
