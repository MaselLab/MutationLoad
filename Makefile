CC = gcc

INCLUDES = -I/home/ulises/gsl/include -I/home/ulises/tskit/c -I/home/ulises/tskit/c/subprojects/kastore/

CFLAGS = $(INCLUDES) -O3 -frounding-math

mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o /home/ulises/tskit/c/tskit/tables.o /home/ulises/tskit/c/subprojects/kastore/kastore.o /home/ulises/tskit/c/tskit/core.o
	$(CC) $(CFLAGS) -O3 -frounding-math pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o tables.o kastore.o core.o -lm -lgsl -lgslcblas -o mutationload

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

dependencies/pcg_basic.o: dependencies/pcg_basic.c
	$(CC) $(CFLAGS) -c dependencies/pcg_basic.c

sharedfunc_flag.o: sharedfunc_flag.c
	$(CC) $(CFLAGS) -c sharedfunc_flag.c

relative_functions.o: relative_functions.c
	$(CC) $(CFLAGS) -c relative_functions.c

absolute_functions.o: absolute_functions.c
	$(CC) $(CFLAGS) -c absolute_functions.c

/home/ulises/tskit/c/tskit/tables.o: /home/ulises/tskit/c/tskit/tables.c
	$(CC) $(CFLAGS) -c -std=c99 /home/ulises/tskit/c/tskit/tables.c

/home/ulises/tskit/c/subprojects/kastore/kastore.o: /home/ulises/tskit/c/subprojects/kastore/kastore.c
	$(CC) $(CFLAGS) -c -std=c99 /home/ulises/tskit/c/subprojects/kastore/kastore.c

/home/ulises/tskit/c/tskit/core.o: /home/ulises/tskit/c/tskit/core.c
	$(CC) $(CFLAGS) -c -std=c99 /home/ulises/tskit/c/tskit/core.c


