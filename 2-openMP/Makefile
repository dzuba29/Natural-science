COMPILER = g++

FLAGS = -std=c++17 -fopenmp
LIBS = -lstdc++fs

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
	$(COMPILER) $(FLAGS) main.cpp -o main $(LIBS)

clean:
	$(RM) *.o main