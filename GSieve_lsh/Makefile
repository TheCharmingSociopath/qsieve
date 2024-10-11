CC=g++
OPTFLAGS= -march=native -O2 -pipe
STRICT= -std=c++17 -pedantic -Wall -fno-builtin
NTL_INCLUDE_PATH=-I/Users/adityamorolia/Projects/pqc-paper/ntl-11.5.1/build/include
NTL_LIB_PATH=-L/Users/adityamorolia/Projects/pqc-paper/ntl-11.5.1/build/lib
GMP_LIB_PATH=-L/opt/homebrew/opt/gmp/lib

## Adi
CFLAGS= $(STRICT) $(OPTFLAGS) NTL_INCLUDE_PATH
LIBS = NTL_LIB_PATH GMP_LIB_PATH -lm -lgmp -lntl
## 

OBJ_DIR = obj
BIN_DIR = bin
BASE = common sampler gsieve
MAIN  = gsieve generate

##
TEMP_BASE_OBJ = $(addsuffix .o, $(BASE)  )
BASE_OBJ = $(addprefix $(OBJ_DIR)/, $(TEMP_BASE_OBJ) )

all: $(MAIN)

$(OBJ_DIR)/%.o: %.cc
	$(CC) $(CONFOPT) $(CFLAGS) -c   $< -o $@

$(MAIN): $(BASE_OBJ)
	$(CC) $@_main.cc $(BASE_OBJ) -o $(BIN_DIR)/$@  $(LIBS) $(CFLAGS)

clean:
	rm bin/* obj/*
