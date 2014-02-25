SRC_DIR= .
CFLAGS+=-w -g -fPIC 
CFLAGS+=$(addprefix -I, $(SRC_DIR))
CPPFILES=$(shell find $(SRC_DIR) -maxdepth 1 -name "*.cpp")
OBJS=$(CPPFILES:%.cpp=%.o)
TARGET=libplot.a

all: $(TARGET)

$(TARGET):$(OBJS)
	ar rcs $(TARGET) $(OBJS) 

$(OBJS):%.o:%.cpp
	g++  -c $(CFLAGS) -o $@ $<

clean:
	rm -f $(addsuffix /*.d, $(SRC_DIR)) $(addsuffix /*.o, $(SRC_DIR)) $(TARGET)
