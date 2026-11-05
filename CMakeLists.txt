cmake_minimum_required(VERSION 3.10)
project(highway_ukf_tracking)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PCL REQUIRED COMPONENTS common io visualization)

add_executable(highway_ukf
  src/main.cpp
  src/ukf.cpp
  src/sim/sensors.cpp
  src/render/render.cpp)
target_include_directories(highway_ukf PRIVATE src ${PCL_INCLUDE_DIRS})
target_compile_definitions(highway_ukf PRIVATE ${PCL_DEFINITIONS}
  PCD_DIR="${CMAKE_SOURCE_DIR}/data/pcd")
target_compile_options(highway_ukf PRIVATE -Wall)
target_link_libraries(highway_ukf PRIVATE Eigen3::Eigen ${PCL_LIBRARIES})
