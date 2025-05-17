LIBOMP := $(shell find /usr/lib/llvm-* -name "libomp.so" | sed 's/libomp.so//')
ifndef LIBOMP
$(error LIBOMP is not set, you need to install libomp-dev)
endif

CXX = g++
AR = ar
CXXFLAGS := -std=c++17 -Wall -pthread -fopenmp
DEBUGFLAG := -g
CPPFLAGS ?= -MMD -MP -mavx2
LDFLAGS := -lpthread -lgmp -lgoldilocks -lssl -lcrypto
ASFLAGS := -felf64

NVCC := /usr/local/cuda/bin/nvcc

INC := -I./include -I./goldilocks/src
LINK := -L./goldilocks/
SRC_DIR := ./src
SRCS := $(shell find $(SRC_DIR) -name '*.cpp')

# g++ -g -std=c++17 -o test test.cpp ./src/* -lssl -lcrypto -lpthread -lgoldilocks -I./include -L./goldilocks/ -fopenmp -mavx2 -I./goldilocks/src -lgmp
LIB_NAME := libgoldilocks.a
LIB_DIR := ./goldilocks
LIB_SRC_DIR := ./goldilocks/src
LIB_SRCS := $(shell find $(LIB_SRC_DIR) -name '*.cpp')
LIB_OBJS := $(LIB_SRCS:%.cpp=%.o)

lib: $(LIB_OBJS)
	$(AR) rcs $(LIB_DIR)/$(LIB_NAME) $(LIB_OBJS)
	$(RM) $(LIB_SRC_DIR)/*.o
	$(RM) $(LIB_SRC_DIR)/*.d

%.cpp.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

test: $(LIB_DIR)/$(LIB_NAME)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o test ./test.cpp $(SRCS) $(LDFLAGS) $(INC) $(LINK) 

clean:
	$(RM) $(LIB_DIR)/$(LIB_NAME)
	$(RM) ./test