TARGET=mpi-1

N_PROC=4
N=8

default: $(TARGET) run

$(TARGET): clean
	mpic++ main.cpp -o $@

clean:
	rm -f $(TARGET)

run:
	mpirun -np $(N_PROC) $(TARGET) 

#-n $(N)
