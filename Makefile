CC = gcc
CFLAGS = -O3 -frounding-math -I/home/mawass/tskit/c
LDFLAGS = -L/home/mawass/tskit/c
LDLIBS = -ltskit -lm -lgsl -lgslcblas

OBJS = main.o \
       dependencies/pcg_basic.o \
       sharedfunc_flag.o \
       relative_functions.o \
       absolute_functions.o

mutationload: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS) mutationload
