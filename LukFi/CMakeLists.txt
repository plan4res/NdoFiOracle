add_executable(LukFiOrcl LukFiOrcl.C test/test.cpp)
target_include_directories(LukFiOrcl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(LukFiOrcl PRIVATE NDOFiOracle::NDOFiOracle CPLEX::Cplex)
