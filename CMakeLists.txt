cmake_minimum_required(VERSION 3.5)
project(g_map)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(GROMACS 5.1 REQUIRED)

add_definitions(${GROMACS_DEFINITIONS})
message(${GROMACS_INCLUDE_DIRS})
include_directories(${GROMACS_INCLUDE_DIRS})

set(SOURCE_FILES src/g_map.cpp)

include_directories(src)

add_executable(g_map ${SOURCE_FILES})
target_link_libraries(g_map ${GROMACS_LIBRARIES})
install(FILES CMakeLists.txt.g_map
        DESTINATION ${DATA_INSTALL_DIR}/g_map
        RENAME CMakeLists.txt
        COMPONENT development)
install(FILES README g_map.cpp Makefile.pkg
        DESTINATION ${DATA_INSTALL_DIR}/g_map
        COMPONENT development)

install(FILES cmake/FindGROMACS.cmake
        DESTINATION ${DATA_INSTALL_DIR}/g_map/cmake
        COMPONENT development)