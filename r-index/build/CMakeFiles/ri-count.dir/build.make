# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/niklas/code/k_index_github/r-index

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/niklas/code/k_index_github/r-index/build

# Include any dependencies generated for this target.
include CMakeFiles/ri-count.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ri-count.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ri-count.dir/flags.make

CMakeFiles/ri-count.dir/ri-count.cpp.o: CMakeFiles/ri-count.dir/flags.make
CMakeFiles/ri-count.dir/ri-count.cpp.o: ../ri-count.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/niklas/code/k_index_github/r-index/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ri-count.dir/ri-count.cpp.o"
	/usr/bin/g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ri-count.dir/ri-count.cpp.o -c /home/niklas/code/k_index_github/r-index/ri-count.cpp

CMakeFiles/ri-count.dir/ri-count.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ri-count.dir/ri-count.cpp.i"
	/usr/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/niklas/code/k_index_github/r-index/ri-count.cpp > CMakeFiles/ri-count.dir/ri-count.cpp.i

CMakeFiles/ri-count.dir/ri-count.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ri-count.dir/ri-count.cpp.s"
	/usr/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/niklas/code/k_index_github/r-index/ri-count.cpp -o CMakeFiles/ri-count.dir/ri-count.cpp.s

CMakeFiles/ri-count.dir/ri-count.cpp.o.requires:

.PHONY : CMakeFiles/ri-count.dir/ri-count.cpp.o.requires

CMakeFiles/ri-count.dir/ri-count.cpp.o.provides: CMakeFiles/ri-count.dir/ri-count.cpp.o.requires
	$(MAKE) -f CMakeFiles/ri-count.dir/build.make CMakeFiles/ri-count.dir/ri-count.cpp.o.provides.build
.PHONY : CMakeFiles/ri-count.dir/ri-count.cpp.o.provides

CMakeFiles/ri-count.dir/ri-count.cpp.o.provides.build: CMakeFiles/ri-count.dir/ri-count.cpp.o


# Object files for target ri-count
ri__count_OBJECTS = \
"CMakeFiles/ri-count.dir/ri-count.cpp.o"

# External object files for target ri-count
ri__count_EXTERNAL_OBJECTS =

ri-count: CMakeFiles/ri-count.dir/ri-count.cpp.o
ri-count: CMakeFiles/ri-count.dir/build.make
ri-count: CMakeFiles/ri-count.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/niklas/code/k_index_github/r-index/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ri-count"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ri-count.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ri-count.dir/build: ri-count

.PHONY : CMakeFiles/ri-count.dir/build

CMakeFiles/ri-count.dir/requires: CMakeFiles/ri-count.dir/ri-count.cpp.o.requires

.PHONY : CMakeFiles/ri-count.dir/requires

CMakeFiles/ri-count.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ri-count.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ri-count.dir/clean

CMakeFiles/ri-count.dir/depend:
	cd /home/niklas/code/k_index_github/r-index/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/niklas/code/k_index_github/r-index /home/niklas/code/k_index_github/r-index /home/niklas/code/k_index_github/r-index/build /home/niklas/code/k_index_github/r-index/build /home/niklas/code/k_index_github/r-index/build/CMakeFiles/ri-count.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ri-count.dir/depend

