find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

include(${PROJECT_SOURCE_DIR}/cmake/add_gtest.cmake)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/matrix)
