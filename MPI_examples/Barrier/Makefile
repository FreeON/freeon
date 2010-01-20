all : master slave

master : master.c
	/usr/local/openmpi-1.3.4-gcc-4.4.2/bin/mpicc -o $@ $^

slave : slave.c
	/usr/local/openmpi-1.3.4-gcc-4.4.2/bin/mpicc -o $@ $^
