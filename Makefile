CC = gcc
INCLUDES = -I/home/gsl/include
CFLAGS = $(INCLUDES) -O3 -frounding-math
LDFLAGS = -L/home/gsl/lib
mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o
	$(CC) $(LDFLAGS) -O3 -frounding-math dependencies/pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o -lm -lgsl -lgslcblas -o mutationload
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
/home/tskit/c/tskit/tables.o: /home/tskit/c/tskit/tables.c
	$(CC) $(CFLAGS) -c -std=c99 /home/tskit/c/tskit/tables.c
/home/tskit/c/subprojects/kastore/kastore.o: /home/tskit/c/subprojects/kastore/kastore.c
	$(CC) $(CFLAGS) -c -std=c99 /home/tskit/c/subprojects/kastore/kastore.c
/home/tskit/c/tskit/core.o: /home/tskit/c/tskit/core.c
	$(CC) $(CFLAGS) -c -std=c99 /home/tskit/c/tskit/core.c
/home/tskit/c/tskit/trees.o: /home/tskit/c/tskit/trees.c
	$(CC) $(CFLAGS) -c -std=c99 /home/tskit/c/tskit/trees.c