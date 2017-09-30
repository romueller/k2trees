cmake_minimum_required(VERSION 3.0)
project(k2trees)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -O3 -msse4.1 -lsdsl")

set(SOURCE_FILES Utility.hpp Utility.cpp K2Tree.hpp RowTree.hpp StaticBasicTree.hpp StaticHybridTree.hpp StaticBasicRectangularTree.hpp StaticUnevenRectangularTree.hpp StaticRowTree.hpp StaticHybridRowTree.hpp StaticMiniRowTree.hpp StaticMiniK2Tree.hpp StaticUnevenRectangularOrMiniTree.hpp)

add_library(k2trees ${SOURCE_FILES})


add_executable(k2trees-rel main.cpp ${SOURCE_FILES})
add_executable(k2trees-deb main.cpp ${SOURCE_FILES})

target_link_libraries(k2trees-rel sdsl)
target_link_libraries(k2trees-deb sdsl)