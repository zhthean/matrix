function(add_gtest)
  set(options OPTIONAL)
  set(oneValueArgs NAME LIBRARY)
  set(multiValueArgs FILES ADDITIONAL_HEADER ADDITIONAL_LIBRARIES )
  cmake_parse_arguments(PARSE_ARGV 0 arg
    "${options}" "${oneValueArgs}" "${multiValueArgs}"
  )

  if(NOT arg_NAME OR NOT arg_FILES)
    message(FATAL_ERROR "NAME and FILES must be specified for add_gtest.")
  endif()

   add_executable(${arg_NAME} ${arg_FILES})

   target_include_directories(
       ${arg_NAME}
       PRIVATE ${arg_ADDITIONAL_HEADER}
   )

   target_link_libraries(
       ${arg_NAME}
       PRIVATE
       ${arg_LIBRARY}
       ${arg_ADDITIONAL_LIBRARIES}
       GTest::gtest
       GTest::gtest_main
   )

   add_test(
       NAME ${arg_NAME}
       COMMAND ${arg_NAME}
   )
endfunction()
