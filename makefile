# Matthew Pautz
# CPSC 1020 001, Sp21
# mpautz@clemson.edu
# Cathy Kittelstad
#
# makefile for programming assignment 2

CC=g++
CFLAGS= -Wall -g -O9
TARGET=probe
FILE=

C_SRCS :=probe.cpp 
HDRS :=  
OBJS :=probe.o 

add: build
	@echo "Done."

build: $(OBJS)
	$(CC) $(OBJS) -o $(TARGET)

# special build rule
%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

run:
	./$(TARGET) $(FILE)

leakcheck: 
	valgrind --leak-check=yes --track-origins=yes ./$(TARGET) $(FILE)

clean:
	rm $(TARGET) $(OBJS)

debug:
	gdb -tui ./$(TARGET) $(FILE)

time:
	time ./$(TARGET)