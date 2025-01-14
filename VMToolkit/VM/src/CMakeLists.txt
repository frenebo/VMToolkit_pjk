

# Find all source and header files in the current directory and its subdirectories
file(GLOB_RECURSE _VMTOOLKIT_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file(GLOB_RECURSE _VMTOOLKIT_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Set these files as accessible in the parent scope
SET(VMTOOLKIT_HDRS ${_VMTOOLKIT_HDRS} PARENT_SCOPE)
SET(VMTOOLKIT_SRCS ${_VMTOOLKIT_SRCS} PARENT_SCOPE)




