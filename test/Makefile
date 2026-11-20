# Host unit tests for BlueIMU (doctest).
#   make          build and run both variants (libm atan2f and BLUEIMU_FAST_TRIG)
#   make exact    only the libm variant
#   make fast     only the fast-trig variant
#   make clean
# Extra doctest options: make ARGS="-s -tc='*gimbal*'"

CXX      ?= g++
STD      ?= -std=gnu++11
LIB_DIR  := ..
BUILD    := build
INCLUDES := -Istub -I$(LIB_DIR) -I.

# The library itself is held to strict warnings, including no implicit double math.
LIB_WARN  := -Wall -Wextra -Wshadow -Wdouble-promotion -Wfloat-conversion -Werror
TEST_WARN := -Wall -Wextra

TESTS := test_math.cpp test_imu.cpp test_euler.cpp

.PHONY: all exact fast clean
all: exact fast

exact: $(BUILD)/tests_exact
	$(abspath $<) $(ARGS)

fast: $(BUILD)/tests_fast
	$(abspath $<) $(ARGS)

$(BUILD):
	mkdir -p $@

# doctest's main is slow to compile; build it once
$(BUILD)/main.o: test_main.cpp vendor/doctest.h | $(BUILD)
	$(CXX) $(STD) -O1 -c $< -o $@

$(BUILD)/lib_exact.o: $(LIB_DIR)/BlueIMU.cpp $(LIB_DIR)/BlueIMU.h | $(BUILD)
	$(CXX) $(STD) -O2 $(LIB_WARN) $(INCLUDES) -c $< -o $@

$(BUILD)/lib_fast.o: $(LIB_DIR)/BlueIMU.cpp $(LIB_DIR)/BlueIMU.h | $(BUILD)
	$(CXX) $(STD) -O2 $(LIB_WARN) -DBLUEIMU_FAST_TRIG $(INCLUDES) -c $< -o $@

$(BUILD)/tests_exact: $(BUILD)/main.o $(BUILD)/lib_exact.o $(TESTS) helpers.h $(LIB_DIR)/BlueIMU.h
	$(CXX) $(STD) -O1 $(TEST_WARN) $(INCLUDES) $(BUILD)/main.o $(BUILD)/lib_exact.o $(TESTS) -o $@ -lm

$(BUILD)/tests_fast: $(BUILD)/main.o $(BUILD)/lib_fast.o $(TESTS) helpers.h $(LIB_DIR)/BlueIMU.h
	$(CXX) $(STD) -O1 $(TEST_WARN) $(INCLUDES) -DBLUEIMU_FAST_TRIG $(BUILD)/main.o $(BUILD)/lib_fast.o $(TESTS) -o $@ -lm

clean:
	rm -rf $(BUILD)
