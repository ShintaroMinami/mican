CC	= gcc
CFLAGS  = -O3 -Wall
LIBS    = -lm
VPATH   = src
SRC     = $(shell ls $(VPATH)/*.c)
HEADERS = $(shell ls $(VPATH)/*.h)
OBJS    = $(SRC:.c=.o)

TARGET  = mican

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

$(OBJS): $(HEADERS)

clean:
	rm -f $(OBJS) $(TARGET) $(TARGETLIB)
