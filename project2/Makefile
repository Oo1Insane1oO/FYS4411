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

cleanData:
	rm data/w1/N2/*
	rm data/w1/N6/*
	rm data/w1/N12/*
	rm data/w1/N20/*
	rm data/w01/N2/*
	rm data/w01/N6/*
	rm data/w01/N12/*
	rm data/w01/N20/*
	rm data/w02/N2/*
	rm data/w02/N6/*
	rm data/w02/N12/*
	rm data/w02/N20/*
	rm data/w03/N2/*
	rm data/w03/N6/*
	rm data/w03/N12/*
	rm data/w03/N20/*
	rm data/w04/N2/*
	rm data/w04/N6/*
	rm data/w04/N12/*
	rm data/w04/N20/*
	rm data/w05/N2/*
	rm data/w05/N6/*
	rm data/w05/N12/*
	rm data/w05/N20/*

.PHONY: clean build

$(EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJECTS): $(OBJ_DIR)/%.o : $(SRC)/%.cpp
	mkdir -p $(@D)
	$(CC) -c $(DEPS) -MMD $(CXXFLAGS) $< -o $@
