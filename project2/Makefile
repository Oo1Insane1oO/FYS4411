CC = g++
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

.PHONY: clean build

$(EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJECTS): $(OBJ_DIR)/%.o : $(SRC)/%.cpp
	mkdir -p $(@D)
	$(CC) -c $(DEPS) -MMD $(CXXFLAGS) $< -o $@
