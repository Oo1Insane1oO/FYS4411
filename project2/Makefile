LDFLAGS = -fopenmp -lpthread
CXXFLAGS = -I eigen/ -O3 -std=c++14

ifeq ($(target),debug)
	LDFLAGS := -fopenmp -lpthread
	CXXFLAGS := -I eigen/ -g -std=c++14 -Wall
endif

vmc_main: src/methods.cpp src/basis.cpp src/vmc.cpp src/tests.cpp src/main.cpp
	g++ $(CXXFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm vmc_main
