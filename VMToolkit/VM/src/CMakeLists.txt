
# macro(add_source_dir dir)
# add_subdirectory(${dir})

# # find all source files in this directory
# file(GLOB VMTOOLKIT_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp)
# file(GLOB VMTOOLKIT_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.hpp)

# # add them to the source group
# source_group(${dir} FILES ${VMTOOLKIT_SRCS} ${VMTOOLKIT_HDRS})

# # add them to the main samos list
# list(APPEND _vmtutorial_sources ${VMTOOLKIT_SRCS} ${VMTOOLKIT_HDRS})
# endmacro(add_source_dir)
# #################


# Find all source and header files in the current directory and its subdirectories
file(GLOB_RECURSE _VMTOOLKIT_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file(GLOB_RECURSE _VMTOOLKIT_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Set these files as accessible in the parent scope
SET(VMTOOLKIT_HDRS ${_VMTOOLKIT_HDRS} PARENT_SCOPE)
SET(VMTOOLKIT_SRCS ${_VMTOOLKIT_SRCS} PARENT_SCOPE)




