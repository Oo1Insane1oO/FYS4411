CC = g++
LDFLAGS	= -fopenmp -lpthread
CXXFLAGS = -I eigen/ -O3 -std=c++14

DEBUG ?= OFF
ifeq ($(DEBUG),ON)
	CXXFLAGS := -I eigen/ -g -std=c++14 -Wall -Wextra
endif

EXEC = vmc_main
SRC = src
OBJ_DIR = bin
SOURCES = $(wildcard $(SRC)/*.cpp)
OBJECTS = $(patsubst $(SRC)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

build: $(EXEC)

clean:
	rm -rf $(OBJ_DIR)

.PHONY: clean build

$(EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJECTS): $(OBJ_DIR)/%.o : $(SRC)/%.cpp
	echo Compiling $<
	mkdir -p $(@D)
	$(CC) -c $(CXXFLAGS) -o $@ $<
