COMPILER = mpic++

FLAGS = -std=c++0x -fopenmp 

#-Wall -Wextra -Wpedantic

OS := $(shell uname)
ifeq ($(OS),Darwin)
	COMPILER = g++-8
endif

ifdef SYSTEMROOT
    RM = del /Q
else
    RM = rm -f
endif


default: all

all:
	$(COMPILER) $(FLAGS) main.cpp -o main

clean:
	$(RM) *.o main