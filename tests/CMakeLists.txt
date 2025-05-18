find_package(GTest REQUIRED)

include(${PROJECT_SOURCE_DIR}/cmake/create_gtest.cmake)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/matrix)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vector)
