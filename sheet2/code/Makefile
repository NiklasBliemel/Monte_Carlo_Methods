CC = gcc
CXX = g++
CXXFLAGS = -O3 -std=c++11 -I.

PROGRAM = main
OBJS = main.o Grid.o Ising.o

$(PROGRAM): $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $^
	
clean:
	rm -f $(PROGRAM) $(OBJS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 
	
%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 
	
%.o: %.c $(HDRS)
	$(CC) $(CXXFLAGS) -c $< 