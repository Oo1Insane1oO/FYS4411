CC = g++
LDFLAGS	= -I eigen/ -fopenmp -lpthread
CXXFLAGS = -O3 -std=c++14

DEBUG ?= OFF
ifeq ($(DEBUG),ON)
	CXXFLAGS := -g -std=c++14 -Wall -Wextra
endif

vmc_main: src/*.cpp
	$(CC) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -rf $(OBJ)

.PHONY: clean
