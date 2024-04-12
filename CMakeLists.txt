CMAKE_MINIMUM_REQUIRED(VERSION 3.6)

PROJECT(UMSC)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# SET(CMAKE_BUILD_TYPE Release)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})

# Add particle system library
FILE(GLOB umsc_inc ${UMSC_SOURCE_DIR}/include/*.h ${UMSC_SOURCE_DIR}/vcc/*.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS}
                    ./include ./vcc)

ADD_EXECUTABLE(umsc ${umsc_inc} simple_viewer.cpp)

TARGET_LINK_LIBRARIES(umsc ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} pthread)
