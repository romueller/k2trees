CXX=g++
CXXFLAGS=-std=c++11 -pthread -O3 -msse4.1
LDFLAGS=-lsdsl
ARFLAGS=rc
BUILD=build
OBJ_DIR=$(BUILD)/objects
TARGET=
LIB_TARGET=libk2trees.a
INCLUDE=
SRC=$(wildcard *.cpp)
HEADERS=$(wildcard *.hpp)

OBJECTS=$(SRC:%.cpp=$(OBJ_DIR)/%.o)

INSTALL_PREFIX?=/usr/local

lib: build $(BUILD)/$(LIB_TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(BUILD)/$(LIB_TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	ar $(ARFLAGS) $@ $(OBJECTS)

.PHONY: lib build install clean

build:
	@mkdir -p $(OBJ_DIR)

clean:
	rm -rf build/*


install:
	mkdir -p $(INSTALL_PREFIX)/include/k2trees
	mkdir -p $(INSTALL_PREFIX)/lib
	cp $(HEADERS) $(INSTALL_PREFIX)/include/k2trees/
	cp $(BUILD)/$(LIB_TARGET) $(INSTALL_PREFIX)/lib/

uninstall:
	rm -rf $(INSTALL_PREFIX)/include/k2trees
	rm $(INSTALL_PREFIX)/lib/$(LIB_TARGET)


