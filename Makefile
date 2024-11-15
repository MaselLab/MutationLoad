
ifeq ($(whereisrun), walid)
CC = gcc
INCLUDES = -I/home/wmawass/gsl/include
CFLAGS = $(INCLUDES) -O3 -frounding-math
LDFLAGS = -L/home/wmawass/gsl/lib
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
endif


ifeq ($(whereisrun), ulises)
mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o
	gcc -O3 -frounding-math dependencies/pcg_basic.o main.o sharedfunc_flag.o  relative_functions.o absolute_functions.o -lm -lgsl -lgslcblas -o mutationload

main.o: main.c
	gcc -O3 -frounding-math -c main.c
endif

ifeq ($(whereisrun), HPC)

CC = gcc

INCLUDES = -I/home/u5/walidmawass/tskit/c -I/home/u5/walidmawass/tskit/c/subprojects/kastore/

CFLAGS = $(INCLUDES) -O3 -frounding-math


mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o /home/u5/walidmawass/tskit/c/tskit/tables.o /home/u5/walidmawass/tskit/c/subprojects/kastore/kastore.o /home/u5/walidmawass/tskit/c/tskit/core.o /home/u5/walidmawass/tskit/c/tskit/trees.o
	$(CC) -O3 -frounding-math pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o tables.o kastore.o core.o trees.o -lm -lgsl -lgslcblas -o mutationload
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
/home/u5/walidmawass/tskit/c/tskit/tables.o: /home/u5/walidmawass/tskit/c/tskit/tables.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u5/walidmawass/tskit/c/tskit/tables.c
/home/u5/walidmawass/tskit/c/subprojects/kastore/kastore.o: /home/u5/walidmawass/tskit/c/subprojects/kastore/kastore.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u5/walidmawass/tskit/c/subprojects/kastore/kastore.c
/home/u5/walidmawass/tskit/c/tskit/core.o: /home/u5/walidmawass/tskit/c/tskit/core.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u5/walidmawass/tskit/c/tskit/core.c
/home/u5/walidmawass/tskit/c/tskit/genotypes.o: /home/u5/walidmawass/tskit/c/tskit/genotypes.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u5/walidmawass/tskit/c/tskit/genotypes.c
/home/u5/walidmawass/tskit/c/tskit/trees.o: /home/u5/walidmawass/tskit/c/tskit/trees.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u5/walidmawass/tskit/c/tskit/trees.c
endif
