cmake_minimum_required(VERSION 2.8.3)
project(grid_path_searcher)
add_compile_options(-std=c++14)

# Set compile option for debug
#set(CMAKE_CXXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
SET(CMAKE_BUILD_TYPE Debug)

find_package(catkin REQUIRED COMPONENTS 
    roscpp 
    std_msgs
    nav_msgs
    visualization_msgs 
)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)

set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

catkin_package(
  INCLUDE_DIRS include
)

include_directories(
    include 
    SYSTEM 
    third_party
    ${catkin_INCLUDE_DIRS} 
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -O0 -Wall") # -Wextra -Werror

add_executable( demo_node 
    src/demo_node.cpp
    src/Astar_searcher.cpp
    src/read_only/JPS_utils.cpp
    src/read_only/JPS_searcher.cpp
    )

target_link_libraries(demo_node 
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES} 
)

add_executable ( random_complex 
    src/random_complex_generator.cpp )

target_link_libraries( random_complex
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES} )  
