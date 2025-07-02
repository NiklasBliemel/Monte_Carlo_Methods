CXX = g++
CXXFLAGS = -O3 -std=c++11 -I.

PROGRAMS = main_1 main_2 main_3 main_4 main_5

OBJS_COMMON = Grid.o Ising.o Routines.o

all: $(PROGRAMS)

main_1: main_1.o $(OBJS_COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $^

main_2: main_2.o $(OBJS_COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $^

main_3: main_3.o $(OBJS_COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $^

main_4: main_4.o $(OBJS_COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $^

main_5: main_5.o $(OBJS_COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(PROGRAMS) *.o
