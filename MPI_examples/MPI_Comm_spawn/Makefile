.PHONY : all clean

all : spawn other

clean :
	rm -f spawn other

spawn : main.c
	$(CC) -o $@ $^

other : child.c
	$(CC) -o $@ $^
