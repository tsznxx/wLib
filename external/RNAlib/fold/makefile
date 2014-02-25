SRC_DIR= .
CFLAGS+=-w -g -fPIC 
CFLAGS+=$(addprefix -I, $(SRC_DIR))
CFILES=$(shell find $(SRC_DIR) -maxdepth 1 -name "*.c")
OBJS=$(CFILES:%.c=%.o)
TARGET=libfold.a

all: $(TARGET)

$(TARGET):$(OBJS)
	ar rcs $(TARGET) $(OBJS) 

$(OBJS):%.o:%.c
	gcc  -c $(CFLAGS) -o $@ $<

clean:
	rm -f $(addsuffix /*.d, $(SRC_DIR)) $(addsuffix /*.o, $(SRC_DIR)) $(TARGET)
