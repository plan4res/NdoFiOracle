# --------------------------------------------------------------------------- #
#    Main CMake file for NDOSolver/FiOracle library                           #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.19)
cmake_policy(VERSION 3.12)
project(NDOFiOracle
        VERSION 4.1
        DESCRIPTION "C++ solvers for NonDifferentiable Optimization problems"
        HOMEPAGE_URL https://gitlab.com/frangio68/ndosolver_fioracle_project
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})
set(modNamespace "NDOFiOracle")

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
#include(CTest)

# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# An option is an ON/OFF user-settable cache variable
option(NDOFiOracle_USE_OSI "Whether ${modName} will use Osi or not." ON)
option(NDOFiOracle_USE_CLP "Whether ${modName} will use Clp or not." ON)

option(NDOFiOracle_READ_PATHS "Whether ${modName} will read locations for
                               dependencies or not." ON)

# ----- Requirements -------------------------------------------------------- #
if (NDOFiOracle_READ_PATHS)
    # Read the dependencies default location from a file written in a
    # bash-style syntax since it is shared with the other hand-written Makefiles,
    # so we define a function to convert these macros written as `lhs = rhs` into
    # their equivalent cmake variants, i.e., `set(lhs rhs)`, avoiding comments.
    macro(makefile_macros_to_cmake file)
        FILE(STRINGS ${file} paths)
        foreach (path IN LISTS paths)
            if (path) # ignore empty lines, i.e., \n
                # makes empty the comment lines, i.e., those start with #...
                string(REGEX REPLACE "^[ \t]*[#].*" "" path ${path})
                if (path) # ... and then ignore them

                    string(REGEX MATCH ".*[=]" lhs "${path}") # retrieve the lhs
                    string(REGEX REPLACE "=" "" lhs ${lhs}) # remove =
                    # remove tab or white spaces from head
                    string(STRIP "${lhs}" lhs)

                    # if it's not already defined before in custom makefile-paths...
                    if (NOT DEFINED ${lhs})
                        # then read its rhs...
                        string(REGEX MATCH "[=].*" rhs "${path}") # retrieve the rhs
                        # remove comments starting with # from tails, if any
                        string(REGEX REPLACE "[#].*" "" rhs ${rhs})
                        string(REGEX REPLACE "=" "" rhs ${rhs}) # remove =
                        string(REGEX REPLACE "\"" "" rhs ${rhs}) # remove ", if any
                        # remove tab or white space from head or tail
                        string(STRIP "${rhs}" rhs)

                        # ... and if it is not a default path in Unix...
                        if (NOT (rhs STREQUAL /usr OR
                                 rhs STREQUAL /usr/local))
                            # ... then set the one defined in makefile-default-paths
                            set("${lhs}" ${rhs})
                            message(STATUS "Specify ${lhs}: " ${rhs})
                        endif ()
                    endif ()

                    unset(lhs)
                    unset(rhs)
                endif ()
            endif ()
        endforeach ()

        unset(path)
        unset(paths)
    endmacro ()

    # firstly, set the custom makefile-paths, if any...
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extlib/makefile-paths")
        makefile_macros_to_cmake(extlib/makefile-paths)
    endif ()
    # ... then set the default other ones ignoring those already previously set
    makefile_macros_to_cmake(extlib/makefile-default-paths)

    # Retrieve the system architecture
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        if (UNIX)
            set(ARCH x86-64)
        else () # Windows
            set(ARCH x64)
        endif ()
    else ()
        set(ARCH x86)
    endif ()
endif ()

if (NDOFiOracle_USE_OSI)
    # if Osi_OsiCpx_FOUND and/or Osi_OsiGrb_FOUND have been set to TRUE before
    # outside of here, then avoid `find_package( )` again
    if ((NOT Osi_OsiCpx_FOUND) AND (NOT Osi_OsiGrb_FOUND))
        # OsiCpx and/or OsiGrb needed for OSIMPSolver
        find_package(Osi OPTIONAL_COMPONENTS OsiCpx OsiGrb)
        if ((NOT Osi_OsiCpx_FOUND) AND (NOT Osi_OsiGrb_FOUND))
            set(NDOFiOracle_USE_OSI OFF)
        endif ()
    endif ()
endif ()

if (NDOFiOracle_USE_CLP)
    # if Clp_FOUND has been set to TRUE before outside of here, then avoid
    # `find_package( )` again
    if (NOT Clp_FOUND)
        # OsiClp needed for CutPlane
        find_package(Clp COMPONENTS OsiClp)
        if (NOT Clp_FOUND)
            set(NDOFiOracle_USE_CLP OFF)
        endif ()
    endif ()
endif ()

# Set the WHICH_OSI_MP macro according to these following priorities:
# - use the given WHICH_OSI_MP in input as CMake option, if any
# - use Gurobi, i.e., 1, if available
# - use Cplex, i.e., 2, if available
# - otherwise, use Clp

# if WHICH_OSI_MP is *not* given in input as CMake option...
# (`NOT DEFINED WHICH_OSI_MP` tests `set( )`, not `add_definitions( )`)
if (NOT DEFINED WHICH_OSI_MP)
    # ... then use OsiGrb in NDOFiOracle by default (as in OSIMPSolver.C:60)
    add_definitions(-DWHICH_OSI_MP=2) # use OsiGrb in NDOFiOracle by default
else () # else if WHICH_OSI_MP is given in input as CMake option...
    if (WHICH_OSI_MP MATCHES "^(CPLEX|cplex|CPX|cpx|1)$")
        add_definitions(-DWHICH_OSI_MP=1) # ... then use OsiCpx in NDOFiOracle
    elseif (WHICH_OSI_MP MATCHES "^(GUROBI|gurobi|GRB|grb|2)$")
        add_definitions(-DWHICH_OSI_MP=2) # ... then use OsiGrb in NDOFiOracle
    endif ()
endif ()

if (NDOFiOracle_USE_OSI)
    # if Cplex is found...
    if (Osi_OsiCpx_FOUND)
        # ... and if WHICH_OSI_MP is *not* given in input as CMake option...
        # (`NOT DEFINED WHICH_OSI_MP` tests `set( )`, not `add_definitions( )`)
        if (NOT DEFINED WHICH_OSI_MP)
            # ... then use OsiCpx in NDOFiOracle
            add_definitions(-DWHICH_OSI_MP=1)
        endif ()
    endif ()

    # but, if also Gurobi is found...
    if (Osi_OsiGrb_FOUND)
        # ... and if WHICH_OSI_MP is *not* given in input as CMake option...
        # (`NOT DEFINED WHICH_OSI_MP` tests `set( )`, not `add_definitions( )`)
        if (NOT DEFINED WHICH_OSI_MP)
            # ... then use OsiGrb in NDOFiOracle
            add_definitions(-DWHICH_OSI_MP=2)
        endif ()
    endif ()
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_17)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
# Note: do not GLOB files here.
target_sources(${modName} PRIVATE
               Bundle/Bundle.C
               LukFi/LukFiOrcl.C
               QPPnltMP/BMinQuad.C
               QPPnltMP/MinQuad.C
               QPPnltMP/QPPnltMP.C
               SubGrad/SubGrad.C)

# When using target_include_directories(), PUBLIC means that any targets
# that link to this target also need that include directory.
# Other options are PRIVATE (only affect the current target, not dependencies),
# and INTERFACE (only needed for dependencies).

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
# This means that if a target finds this library in its build directory
# will look into the BUILD_INTERFACE path for its headers, if it finds it
# installed in the system will look into the INSTALL_INTERFACE path.
target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Bundle>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LukFi>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MPTester>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/NDOSlver>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/QPPnltMP>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SubGrad>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Deflection>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Stepsize>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
if (NDOFiOracle_USE_CLP)
    message(STATUS "${modName}: OsiClp found, configuring CutPlane")

    target_sources(${modName} PRIVATE CutPlane/CutPlane.C)
    target_include_directories(
            ${modName} PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CutPlane>)
    target_link_libraries(${modName} PUBLIC Coin::OsiClp)
endif ()

if (NDOFiOracle_USE_OSI)
    target_sources(${modName} PRIVATE OSIMPSolver/OSIMPSolver.C)
    target_include_directories(
            ${modName} PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OSIMPSolver>)

    if (Osi_OsiCpx_FOUND)
        message(STATUS "${modName}: OsiCpx found, configuring OSIMPSolver")
        target_link_libraries(${modName} PUBLIC Coin::OsiCpx)
    endif ()

    if (Osi_OsiGrb_FOUND)
        message(STATUS "${modName}: OsiGrb found, configuring OSIMPSolver")
        target_link_libraries(${modName} PUBLIC Coin::OsiGrb)
    endif ()
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- Subdirectories ------------------------------------------------------ #
#add_subdirectory(TestFi)
#add_subdirectory(LukFi)

if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils/OPTtypes.h
        ${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils/OPTUtils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils/OPTvect.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Bundle/Bundle.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Bundle/MPSolver.h
        ${CMAKE_CURRENT_SOURCE_DIR}/LukFi/LukFiOrcl.h
        ${CMAKE_CURRENT_SOURCE_DIR}/MPTester/MPTester.h
        ${CMAKE_CURRENT_SOURCE_DIR}/NDOSlver/NDOSlver.h
        ${CMAKE_CURRENT_SOURCE_DIR}/NDOSlver/FiOracle.h
        ${CMAKE_CURRENT_SOURCE_DIR}/QPPnltMP/BMinQuad.h
        ${CMAKE_CURRENT_SOURCE_DIR}/QPPnltMP/MinQuad.h
        ${CMAKE_CURRENT_SOURCE_DIR}/QPPnltMP/QPPnltMP.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/SubGrad.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Stepsize.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Deflection.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Deflection/Volume.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Deflection/PrimalDual.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Stepsize/ColorTV.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Stepsize/FumeroTV.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGrad/Stepsize/Polyak.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

if (NDOFiOracle_USE_OSI)
    install(FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/CutPlane/CutPlane.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})
endif ()

if (NDOFiOracle_USE_CLP)
    install(FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/OSIMPSolver/OSIMPSolver.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})
endif ()

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Shipping custom find modules should not be encouraged because third
# party developers should provide for them.
file(GLOB findMods ${CMAKE_CURRENT_LIST_DIR}/cmake/Find*)
install(FILES ${findMods}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the README and LICENSE files.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

# --------------------------------------------------------------------------- #
# Remove from the search paths the cmake directory we added at the beginning.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #
