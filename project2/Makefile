CC = mpic++
LDFLAGS	= -fopenmp -lpthread
CXXFLAGS = -I eigen/ -O3 -fopenmp -std=c++14

DEBUG ?= OFF
ifeq ($(DEBUG),ON)
	CXXFLAGS := -I eigen/ -g -std=c++14 -Wall -Wextra
endif

EXEC = vmc_main
SRC = src
OBJ_DIR = bin
SOURCES = $(wildcard $(SRC)/*.cpp)
OBJECTS = $(patsubst $(SRC)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
DEPS = $($(OBJECTS):.o=.d)

-include $(DEPS)

build: $(EXEC)

clean:
	rm -rf $(OBJ_DIR)
	rm $(EXEC)

createDataDir:
	mkdir data
	mkdir data/w010
	mkdir data/w01
	mkdir data/w05
	mkdir data/w005
	mkdir data/w001
	mkdir data/w010/N2
	mkdir data/w010/N6
	mkdir data/w010/N12
	mkdir data/w010/N20
	mkdir data/w01/N2
	mkdir data/w01/N6
	mkdir data/w01/N12
	mkdir data/w01/N20
	mkdir data/w05/N2
	mkdir data/w05/N6
	mkdir data/w05/N12
	mkdir data/w05/N20
	mkdir data/w005/N2
	mkdir data/w005/N6
	mkdir data/w005/N12
	mkdir data/w005/N20
	mkdir data/w001/N2
	mkdir data/w001/N6
	mkdir data/w001/N12
	mkdir data/w001/N20


cleanData:
	rm -f data/w010/N2/*
	rm -f data/w010/N6/*
	rm -f data/w010/N12/*
	rm -f data/w010/N20/*
	rm -f data/w01/N2/*
	rm -f data/w01/N6/*
	rm -f data/w01/N12/*
	rm -f data/w01/N20/*
	rm -f data/w02/N2/*
	rm -f data/w02/N6/*
	rm -f data/w02/N12/*
	rm -f data/w02/N20/*
	rm -f data/w03/N2/*
	rm -f data/w03/N6/*
	rm -f data/w03/N12/*
	rm -f data/w03/N20/*
	rm -f data/w04/N2/*
	rm -f data/w04/N6/*
	rm -f data/w04/N12/*
	rm -f data/w04/N20/*
	rm -f data/w05/N2/*
	rm -f data/w05/N6/*
	rm -f data/w05/N12/*
	rm -f data/w05/N20/*

.PHONY: clean build

$(EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJECTS): $(OBJ_DIR)/%.o : $(SRC)/%.cpp
	mkdir -p $(@D)
	$(CC) -c $(DEPS) -MMD $(CXXFLAGS) $< -o $@
