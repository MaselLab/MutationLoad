HOME_DIR = /home/u6/micailamarcelle

CC = gcc
INCLUDES = -I$(HOME_DIR)/tskit/c -I$(HOME_DIR)/tskit/c/subprojects/kastore/
CFLAGS = $(INCLUDES) -O3 -frounding-math

mutationload: main.o dependencies/pcg_basic.o sharedfunc_flag.o relative_functions.o absolute_functions.o $(HOME_DIR)/tskit/c/tskit/tables.o $(HOME_DIR)/tskit/c/subprojects/kastore/kastore.o $(HOME_DIR)/tskit/c/tskit/core.o $(HOME_DIR)/tskit/c/tskit/genotypes.o $(HOME_DIR)/tskit/c/tskit/trees.o
	$(CC) -O3 -frounding-math pcg_basic.o main.o sharedfunc_flag.o relative_functions.o absolute_functions.o tables.o kastore.o core.o genotypes.o trees.o -lm -lgsl -lgslcblas -o mutationload
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
$(HOME_DIR)/tskit/c/tskit/tables.o: $(HOME_DIR)/tskit/c/tskit/tables.c
	$(CC) $(CFLAGS) -c -std=c99 $(HOME_DIR)/tskit/c/tskit/tables.c
$(HOME_DIR)/tskit/c/subprojects/kastore/kastore.o: $(HOME_DIR)/tskit/c/subprojects/kastore/kastore.c
	$(CC) $(CFLAGS) -c -std=c99 $(HOME_DIR)/tskit/c/subprojects/kastore/kastore.c
$(HOME_DIR)/tskit/c/tskit/core.o: $(HOME_DIR)/tskit/c/tskit/core.c
	$(CC) $(CFLAGS) -c -std=c99 $(HOME_DIR)/tskit/c/tskit/core.c
$(HOME_DIR)/tskit/c/tskit/genotypes.o: $(HOME_DIR)/tskit/c/tskit/genotypes.c
	$(CC) $(CFLAGS) -c -std=c99 $(HOME_DIR)/tskit/c/tskit/genotypes.c
$(HOME_DIR)/tskit/c/tskit/trees.o: $(HOME_DIR)/tskit/c/tskit/trees.c
	$(CC) $(CFLAGS) -c -std=c99 $(HOME_DIR)/tskit/c/tskit/trees.c

.PHONY: clean
clean:
	rm -f absolute_functions.o 
	rm -f core.o 
	rm -f genotypes.o 
	rm -f kastore.o 
	rm -f main.o 
	rm -f pcg_basic.o 
	rm -f relative_functions.o 
	rm -f shared_func_flag.o 
	rm -f tables.o 
	rm -f trees.o 
	rm -f mutationload
