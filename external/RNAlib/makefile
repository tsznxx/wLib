SRC_DIR= fold plot
CFLAGS+=-w -g -fPIC 
CFLAGS+=$(addprefix -I, $(SRC_DIR))
CFILES=$(shell find $(SRC_DIR) -maxdepth 1 -name "*.c")
COBJS=$(CFILES:%.c=%.o)
CPPFILES=$(shell find $(SRC_DIR) -maxdepth 1 -name "*.cpp")
CPPOBJS=$(CPPFILES:%.cpp=%.o)
TARGET=libRNA.a

all: $(TARGET)

$(TARGET): $(COBJS) $(CPPOBJS)
	ar rcs $(TARGET) $(COBJS) $(CPPOBJS) 

$(COBJS):%.o:%.c
	gcc -c $(CFLAGS) -o $@ $<

$(CPPOBJS):%.o:%.cpp
	g++ -c $(CFLAGS) -o $@ $<

clean:
	rm -f $(addsuffix /*.d, $(SRC_DIR)) $(addsuffix /*.o, $(SRC_DIR)) $(TARGET)
