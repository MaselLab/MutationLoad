CC = gcc

# Add include paths AND the library path to CFLAGS
CFLAGS = -O3 -frounding-math -I/home/mawass/tskit -I/home/mawass/tskit/c -I/home/mawass/tskit/c/subprojects/kastore -L/path/to/your/tskit/lib

# Add the tskit and kastore libraries to LDLIBS
# You will likely need kastore as well, as it's a tskit dependency
LDLIBS = -lm -lgsl -lgslcblas -ltskit -lkastore

OBJS = main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o 
mutationload: $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@
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
/home/mawass/tskit/c/tskit/tables.o: /home/mawass/tskit/c/tskit/tables.c
	$(CC) $(CFLAGS) -c -std=c99 /home/mawass/tskit/c/tskit/tables.c
/home/mawass/tskit/c/subprojects/kastore/kastore.o: /home/mawass/tskit/c/subprojects/kastore/kastore.c
	$(CC) $(CFLAGS) -c -std=c99 /home/mawass/tskit/c/subprojects/kastore/kastore.c
/home/mawass/tskit/c/tskit/core.o: /home/mawass/tskit/c/tskit/core.c
	$(CC) $(CFLAGS) -c -std=c99 /home/mawass/tskit/c/tskit/core.c
/home/mawass/tskit/c/tskit/trees.o: /home/mawass/tskit/c/tskit/trees.c
	$(CC) $(CFLAGS) -c -std=c99 /home/mawass/tskit/c/tskit/trees.c
