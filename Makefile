CC = gcc
FLAGS = -std=c99 -Wall
OPT = -O3
SRC = main.c simplex.c
BIN = test

test: main.c simplex.c
	$(CC) $(FLAGS) $(OPT) -o $(BIN) $(SRC)

clean:
	rm $(BIN)
