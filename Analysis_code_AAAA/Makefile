ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
Target=aaAnom_signal
all:$(Target)
	
aaAnom_signal: aaAnom_signal.cxx
	g++ -o $@ aaAnom_signal.cxx $(ROOTFLAGS) $(ROOTLIBS)
clean:
	rm -f $(Target) 
