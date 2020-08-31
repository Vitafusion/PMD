so:
	gcc -fPIC  -I/usr/share/R/include -c normal_appro.c
	gcc -shared -L/Home/admin1/gsl/lib normal_appro.o -lgsl -lgslcblas -lm -o normal_appro.so
clean:
	rm *.o
	rm *.out
	rm *.so
exp:
	export LD_LIBRARY_PATH
	LD_LIBRARY_PATH=/Home/admin1/gsl/lib
