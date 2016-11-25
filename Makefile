fronhofer_dendritic_metapop : fronhofer_dendritic_metapop.cpp include/classes.h include/procedures.h
	g++ -Wall -g -o fronhofer_dendritic_metapop fronhofer_dendritic_metapop.cpp -lgsl -lgslcblas -lm
	
all:	fronhofer_dendritic_metapop

clean:
	rm -f fronhofer_dendritic_metapop
