CXX := g++
AR := ar
CXXFLAGS := -std=c++17 -g 
LDFLAGS := 
SRC := $(wildcard src/*.cpp)
SRC += ./golilocks/src/*.cpp
OBJ := $(SRC:.cpp=.o)
LIB := 
TEST_SRC := $(wildcard test/test_*.cpp)
TEST_BIN := $(TEST_SRC:.cpp=)

all: $(LIB)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# test: all test.cpp
# 	$(CXX) $(CXXFLAGS) test.cpp -L. -lsumcheck -o test
# 	@echo "\n"
# 	./test

refresh: clean all


test: refresh ${TEST_BIN}

test/%: test/%.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $< -L. -lsumcheck -o $@
	@echo "Built $@"

clean:
	rm -f $(OBJ) $(LIB) $(TEST_BIN)

.PHONY:all