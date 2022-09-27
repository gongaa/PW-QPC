CXX      := g++
# CXXFLAGS := -pedantic-errors -Wall -std=c++17
CXXFLAGS := -std=c++17
LDFLAGS  := -L/usr/lib -lstdc++ -lm
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
TARGET   := program
INCLUDE  := -Iinclude/
IGNORE   := src/Algo/PRNG_MT19937_simd.cpp src/Algo/Event_generator_fast.cpp src/Code/Polar_code.cpp src/Decoder/Decoder_polar.cpp
SRC      :=                        \
   $(wildcard src/Algo/*.cpp) 	   \
   $(wildcard src/Channel/*.cpp)   \
   $(wildcard src/Code/*.cpp) 	   \
   $(wildcard src/Decoder/*.cpp)   \
   $(wildcard src/Encoder/*.cpp)   \
   $(wildcard src/Test/*.cpp)      \
   $(wildcard src/*.cpp)           \
   main.cpp
SRC 	 := $(filter-out $(IGNORE), $(SRC)) 

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES \
         := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
