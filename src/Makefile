CC	= gcc
CFLAGS  = -O3 -Wall
LIBS    = -lm
SRC     = $(shell ls *.c)
HEADERS = $(shell ls *.h)
OBJS    = $(SRC:.c=.o)

TARGET = mican

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

$(OBJS): $(HEADERS)

clean:
	rm -f $(OBJS) $(TARGET)
