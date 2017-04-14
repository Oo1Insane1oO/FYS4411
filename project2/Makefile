CC = g++
LDFLAGS	= -fopenmp -lpthread
CXXFLAGS = -I eigen/ -O3 -std=c++14

DEBUG ?= OFF
ifeq ($(DEBUG),ON)
	CXXFLAGS := -I eigen/ -g -std=c++14 -Wall -Wextra
endif

EXEC = vmc_main
SRC = src

$(EXEC): $(SRC)/*.cpp
	$(CC) $(CXXFLAGS) -o $@  $^ $(LDFLAGS)

clean:
	rm $(EXEC)

.PHONY: clean
