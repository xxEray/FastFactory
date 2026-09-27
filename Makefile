CXX = g++
CXXFLAGS = -O3 -march=native -fopenmp -std=c++17
TARGET = solution

all: $(TARGET)

$(TARGET): src/solution.cpp src/bigint.h src/common.h src/mulmod.h src/speed.h
	$(CXX) $(CXXFLAGS) src/solution.cpp -o $(TARGET)

verify: src/verify.cpp src/bigint.h src/common.h src/mulmod.h src/speed.h
	$(CXX) $(CXXFLAGS) src/verify.cpp -o verify

clean:
	rm -f $(TARGET) verify