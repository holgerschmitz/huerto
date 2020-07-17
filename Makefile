
TARGET=schnar_tests

#OFLAGS  = -g -O0 -Wall
OFLAGS  = -O3 -Wall

INCLUDE = -I/usr/local/include -I/usr/include/hdf5/mpich
CXX     = mpiCC

CXXFLAGS = $(OFLAGS) -std=c++11 -DHUERTO_THREE_DIM

BUILD_DIR = build
BIN_DIR = bin
SOURCES =  \
	io/table_data_source.cpp \
	maths/random.cpp \
	tables/table_lookup.cpp \
	tests/main.cpp \
	tests/maths/test_interpolate1d.cpp \
	tests/maths/test_interpolate2d.cpp \
	tests/maths/vector1d.cpp \
	tests/maths/vector2d.cpp \
	tests/maths/vector3d.cpp \
	tests/io/test_table_data_source.cpp \
	tests/tables/test_table_lookup.cpp

OBJECTS = $(addprefix $(BUILD_DIR)/,$(patsubst %.cpp,%.o,$(SOURCES)))

LDFLAGS = 

LOADLIBS = -lschnek -lm

FULLTARGET = $(BIN_DIR)/$(TARGET)

.PHONY: $(TARGET)
$(TARGET): $(BIN_DIR)/$(TARGET)

# Create parallel directory tree to hold the build files
# from
# http://ismail.badawi.io/blog/2017/03/28/automatic-directory-creation-in-make/
.PRECIOUS: $(BUILD_DIR)/. $(BUILD_DIR)%/.

$(BUILD_DIR)/.:
	mkdir -p $@

$(BUILD_DIR)%/.:
	mkdir -p $@

.SECONDEXPANSION:

$(BUILD_DIR)/%.o: %.cpp | $$(@D)/.
	$(CXX) -o $@ -c $(CXXFLAGS) $(INCLUDE) $<


 
all: $(FULLTARGET)

$(FULLTARGET): $(OBJECTS) 
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)


%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f $(OBJECTS) core $(FULLTARGET)


