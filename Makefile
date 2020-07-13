Make:
	g++ -std=c++11 isochrones.cpp -o isoch -O3 -L /usr/lib/x86_64-linux-gnu/ -lboost_graph -lboost_program_options
	./isoch --help
