ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
FASTJET=`/home/cristian/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
Target=aaAnom_signal
all:$(Target)
	
aaAnom_signal: aaAnom_signal.cxx
	g++ -o $@ aaAnom_signal.cxx $(ROOTFLAGS) $(ROOTLIBS) $(FASTJET)
clean:
	rm -f $(Target) 
