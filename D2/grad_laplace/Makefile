#CC=g++
CC=clang++
CFLAGS=-I.
DEPS=$(wildcard *.h)
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.c=.o)
TARGET=main_laplace.x


all: info $(TARGET)
					$(info You are ready to execute!!)

.PHONY: info

info:
	$(info --------------------------- deps)
	$(info $(DEPS))
	$(info ------------------------ srcs)
	$(info $(SRCS))
	$(info --------------------------- OBJS)
	$(info $(OBJS))

%.o: %.c $(DEPS)
	echo "$<" # frist element of r.h.s
	echo "$^" # all r.h.s.
	echo "$@" # target name (l.h.s.)
	$(CC) -c $< $(CFLAGS) -o $@


$(TARGET): $(OBJS)
	@echo "$^"
	$(CC) -o $@ $(OBJS)


clean:
	$(info ********** cleaning up *********)
	@rm -f *~ *.o $(TARGET) > /dev/null 2>&1
	#@rm -f ../data/*.dat


.PHONY: clean
# try to modify a .c file and recompile the library
