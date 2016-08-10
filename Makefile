all:	inflight.cxx sterile_flux.h sterile_flux.cxx fourmomentum.h fourmomentum.cxx channel.h channel.cxx
	g++ -g -std=c++11  -c inflight.cxx -o inflight.o -I. 
	g++ -g -std=c++11  -c fourmomentum.cxx -o fourmomentum.o -I.
	g++ -g -std=c++11  -c sterile_flux.cxx -o sterile_flux.o -I.  
	g++ -g -std=c++11  -c channel.cxx -o channel.o -I. 
	g++ -g -std=c++11  -c detector.cxx -o detector.o -I.  
	g++ -g -std=c++11  -o inflight inflight.o detector.o fourmomentum.o sterile_flux.o channel.o -lgsl -lgslcblas
	rm *.o

