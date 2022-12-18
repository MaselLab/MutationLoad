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


ifeq ($(whereisrun), uhmclocal)
mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o tables.o kastore.o core.o
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math dependencies/pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o tables.o kastore.o core.o -lm -lgsl -lgslcblas -o mutationload

main.o: main.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c main.c

dependencies/pcg_basic.o: dependencies/pcg_basic.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c dependencies/pcg_basic.c

sharedfunc_flag.o: sharedfunc_flag.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c sharedfunc_flag.c

relative_functions.o: relative_functions.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c relative_functions.c

absolute_functions.o: absolute_functions.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c absolute_functions.c

tables.o: ../tskit/c/tskit/tables.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c -std=c99 ../tskit/c/tskit/tables.c

kastore.o: ../tskit/c/subprojects/kastore/kastore.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c -std=c99 ../tskit/c/subprojects/kastore/kastore.c

core.o: ../tskit/c/tskit/core.c
	gcc -I../tskit/c -I../tskit/c/subprojects/kastore/ -O3 -frounding-math -c -std=c99 ../tskit/c/tskit/core.c

endif


ifeq ($(whereisrun), uhmcfusion)
CC = gcc
INCLUDES = -I/home/uhernandez/gsl/include
CFLAGS = $(INCLUDES) -O3 -frounding-math
LDFLAGS = -L/home/uhernandez/gsl/lib
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


ifeq ($(whereisrun), uhmchpc)
CC = gcc
INCLUDES = -I/home/u12/uliseshmc/tskit/c -I/home/u12/uliseshmc/tskit/c/subprojects/kastore/
CFLAGS = $(INCLUDES) -O3 -frounding-math
mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o /home/u12/uliseshmc/tskit/c/tskit/tables.o /home/u12/uliseshmc/tskit/c/subprojects/kastore/kastore.o /home/u12/uliseshmc/tskit/c/tskit/core.o
	$(CC) -O3 -frounding-math dependencies/pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o tables.o kastore.o core.o -lm -lgsl -lgslcblas -o mutationload
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
/home/u12/uliseshmc/tskit/c/tskit/tables.o: /home/u12/uliseshmc/tskit/c/tskit/tables.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u12/uliseshmc/tskit/c/tskit/tables.c
/home/u12/uliseshmc/tskit/c/subprojects/kastore/kastore.o: /home/u12/uliseshmc/tskit/c/subprojects/kastore/kastore.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u12/uliseshmc/tskit/c/subprojects/kastore/kastore.c
/home/u12/uliseshmc/tskit/c/tskit/core.o: /home/u12/uliseshmc/tskit/c/tskit/core.c
	$(CC) $(CFLAGS) -c -std=c99 /home/u12/uliseshmc/tskit/c/tskit/core.c
endif

ifeq ($(whereisrun), clean)
clean:
	rm ./mutationload -f
	rm *.o -f
	rm -r datafor*
endif
