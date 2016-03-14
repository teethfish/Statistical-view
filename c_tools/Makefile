## DEPENDENCIES ##
HDF5_DIR = /opt/hdf5
CGNS_DIR = /opt/cgns

## COMPILERS
GCC = gcc

SIM_DIR = sim
SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

SRCC =	main.c 		\
	cgns_reader.c	\
        dataproc_3d.c	\
	dataproc_2d.c	\
	tool.c
        
EXTRA = Makefile	\
	main.h		\
	cgns_reader.h	\
	dataproc_3d.h	\
	dataproc_2d.h	\
	tool.h

# compile normally:
# make -g -O2 to optimize
all: COPT += -O2
all: main

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: main

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

main: $(OBJS) 
	$(GCC) $(COPT) $+ $(LDLIBS) -lstdc++ -o $(SIM_DIR)/main

clean:
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/main
