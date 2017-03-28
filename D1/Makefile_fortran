test: test.o random_gen.o Makefile
	gfortran -o test  random_gen.o test.o -lstdc++ -O0 -ggdb3

random_gen.o: random_gen.cpp Makefile
	g++ -o random_gen.o -c random_gen.cpp  -Wall -Wall -O0 -ggdb3

test.o: test.f Makefile
	gfortran -c -o test.o test.f -O0 -ggdb3

clean:
	rm test.o random_gen.o test
