# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build

# Include any dependencies generated for this target.
include CMakeFiles/option_pricing.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/option_pricing.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/option_pricing.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/option_pricing.dir/flags.make

CMakeFiles/option_pricing.dir/wrappers.cpp.o: CMakeFiles/option_pricing.dir/flags.make
CMakeFiles/option_pricing.dir/wrappers.cpp.o: ../wrappers.cpp
CMakeFiles/option_pricing.dir/wrappers.cpp.o: CMakeFiles/option_pricing.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/option_pricing.dir/wrappers.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/option_pricing.dir/wrappers.cpp.o -MF CMakeFiles/option_pricing.dir/wrappers.cpp.o.d -o CMakeFiles/option_pricing.dir/wrappers.cpp.o -c /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/wrappers.cpp

CMakeFiles/option_pricing.dir/wrappers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/option_pricing.dir/wrappers.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/wrappers.cpp > CMakeFiles/option_pricing.dir/wrappers.cpp.i

CMakeFiles/option_pricing.dir/wrappers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/option_pricing.dir/wrappers.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/wrappers.cpp -o CMakeFiles/option_pricing.dir/wrappers.cpp.s

# Object files for target option_pricing
option_pricing_OBJECTS = \
"CMakeFiles/option_pricing.dir/wrappers.cpp.o"

# External object files for target option_pricing
option_pricing_EXTERNAL_OBJECTS =

option_pricing.cpython-310-x86_64-linux-gnu.so: CMakeFiles/option_pricing.dir/wrappers.cpp.o
option_pricing.cpython-310-x86_64-linux-gnu.so: CMakeFiles/option_pricing.dir/build.make
option_pricing.cpython-310-x86_64-linux-gnu.so: CMakeFiles/option_pricing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module option_pricing.cpython-310-x86_64-linux-gnu.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/option_pricing.dir/link.txt --verbose=$(VERBOSE)
	/usr/bin/strip /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build/option_pricing.cpython-310-x86_64-linux-gnu.so

# Rule to build all files generated by this target.
CMakeFiles/option_pricing.dir/build: option_pricing.cpython-310-x86_64-linux-gnu.so
.PHONY : CMakeFiles/option_pricing.dir/build

CMakeFiles/option_pricing.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/option_pricing.dir/cmake_clean.cmake
.PHONY : CMakeFiles/option_pricing.dir/clean

CMakeFiles/option_pricing.dir/depend:
	cd /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11 /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11 /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build /mnt/c/Users/luiso/OneDrive/Documents/github_documents/option_pricing_pybind11/build/CMakeFiles/option_pricing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/option_pricing.dir/depend

