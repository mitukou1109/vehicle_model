cmake_minimum_required(VERSION 3.8)
project(vehicle_model)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(semitrailer_model SHARED
  src/semitrailer_model.cpp
)
target_include_directories(semitrailer_model
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)

install(
  TARGETS semitrailer_model
  EXPORT ${PROJECT_NAME}-export
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION ${EIGEN3_INCLUDE_DIR}
)

install(
  DIRECTORY include/
  DESTINATION include
)

install(
  EXPORT ${PROJECT_NAME}-export
  FILE ${PROJECT_NAME}-config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION share/${PROJECT_NAME}/cmake/
)
