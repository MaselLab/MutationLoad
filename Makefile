CC = gcc

INCLUDES = -I/home/wmawass/gsl/include -I/home/wmawass/tskit/c -I/home/wmawass/tskit/c/subprojects/kastore/

CFLAGS = $(INCLUDES) -O3 -frounding-math

LDFLAGS = -L/home/wmawass/gsl/lib

mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o /home/wmawass/tskit/c/tskit/tables.o /home/wmawass/tskit/c/subprojects/kastore/kastore.o /home/wmawass/tskit/c/tskit/core.o
	$(CC) $(LDFLAGS) -O3 -frounding-math pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o tables.o kastore.o core.o -lm -lgsl -lgslcblas -o mutationload

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

/home/wmawass/tskit/c/tskit/tables.o: /home/wmawass/tskit/c/tskit/tables.c
	$(CC) $(CFLAGS) -c -std=c99 /home/wmawass/tskit/c/tskit/tables.c

/home/wmawass/tskit/c/subprojects/kastore/kastore.o: /home/wmawass/tskit/c/subprojects/kastore/kastore.c
	$(CC) $(CFLAGS) -c -std=c99 /home/wmawass/tskit/c/subprojects/kastore/kastore.c

/home/wmawass/tskit/c/tskit/core.o: /home/wmawass/tskit/c/tskit/core.c
	$(CC) $(CFLAGS) -c -std=c99 /home/wmawass/tskit/c/tskit/core.c


