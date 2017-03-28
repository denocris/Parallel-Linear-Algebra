test: test.o inverse_laplace_operator.o Makefile
	gfortran -o test  inverse_laplace_operator.o test.o -lstdc++ -O0 -ggdb3

inverse_laplace_operator.o: inverse_laplace_operator.cpp Makefile
	g++ -o inverse_laplace_operator.o -c inverse_laplace_operator.cpp  -Wall -Wall -O0 -ggdb3

test.o: test.f Makefile
	gfortran -c -o test.o test.f -O0 -ggdb3

clean:
	rm test.o inverse_laplace_operator.o test
