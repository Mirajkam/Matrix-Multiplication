CC=g++
CFLAGS= -O2 -Wall -fopenmp 
LDFLAGS=

#make the program
hw2: hw2.cpp
	$(CC) $(CFLAGS) -o $@ $?

hw3: hw3.cpp
	$(CC) $(CFLAGS) -o $@ $?

#cleanup function
clean1:
	rm hw2

clean2:
	rm hw3
	

#-gencode arch=compute_20,code=sm_20
