cmake_minimum_required (VERSION 3.12)

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

option(ENABLE_TESTS "Enable or disable test cases" OFF)

project ("matrix")

include_directories(${PROJECT_SOURCE_DIR}/include)

# Add source to this project's executable.
add_executable (matrix "src/main.cpp" "src/matrix.cpp"  "src/utils/matrix_func.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET matrix PROPERTY CXX_STANDARD 20)
endif()


if (ENABLE_TESTS)
  enable_testing()
  add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
endif()
