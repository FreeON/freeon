.PHONY : all clean

all : spawn other

clean :
	rm -f spawn other

spawn : main.c
	$(CC) -o $@ -O0 -g $^

other : child.c
	$(CC) -o $@ -O0 -g $^
