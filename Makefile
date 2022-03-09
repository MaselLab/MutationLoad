mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o
	gcc -L/home/wmawass/gsl/lib -O3 -frounding-math dependencies/pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o -lm -lgsl -lgslcblas -o mutationload

main.o: main.c
	gcc -I/home/wmawass/gsl/include -O3 -frounding-math -c main.c

dependencies/pcg_basic.o: dependencies/pcg_basic.c
	gcc -I/home/wmawass/gsl/include -O3 -frounding-math -c dependencies/pcg_basic.c

sharedfunc_flag.o: sharedfunc_flag.c
	gcc -I/home/wmawass/gsl/include -O3 -frounding-math -c sharedfunc_flag.c

relative_functions.o: relative_functions.c
	gcc -I/home/wmawass/gsl/include -O3 -frounding-math -c relative_functions.c

absolute_functions.o: absolute_functions.c
	gcc -I/home/wmawass/gsl/include -O3 -frounding-math -c absolute_functions.c

