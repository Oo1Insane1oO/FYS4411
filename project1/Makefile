LDFLAGS = -fopenmp -lpthread
CXXFLAGS = -I eigen/ -O3 -std=c++11 -Wall

main: src/Coulomb_Functions.cpp src/methods.cpp src/basis.cpp src/tests.cpp src/main.cpp
	g++ $(CXXFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm main
