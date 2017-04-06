LDFLAGS = -fopenmp -lpthread
CXXFLAGS = -I eigen/ -O3 -std=c++11 -Wall

vmc_main: src/methods.cpp src/basis.cpp src/vmc.cpp src/tests.cpp src/main.cpp
	g++ $(CXXFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm main
