add_executable(TestFi TestFi.C test/test.cpp)
target_include_directories(TestFi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TestFi PRIVATE NDOFiOracle::NDOFiOracle)
