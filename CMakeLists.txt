cmake_minimum_required(VERSION 3.8)
project(vehicle_model)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake_auto REQUIRED)

find_package(Eigen3 REQUIRED NO_MODULE)

ament_auto_find_build_dependencies()

ament_auto_add_library(semitrailer_model SHARED
  src/semitrailer_model.cpp
)
target_include_directories(semitrailer_model
  PUBLIC
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
)

ament_export_include_directories(
  ${EIGEN3_INCLUDE_DIR}
)

ament_auto_package()
