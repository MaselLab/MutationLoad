mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o
	gcc -O3 -frounding-math dependencies/pcg_basic.o main.o sharedfunc_flag.o  relative_functions.o -lm -lgsl -lgslcblas -o mutationload

main.o: main.c
	gcc -O3 -frounding-math -c main.c

