cmake_minimum_required(VERSION 3.18)

project(MuSpinDynamic VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(MUSPINDYNAMIC_IS_TOP_LEVEL ON)
else()
    set(MUSPINDYNAMIC_IS_TOP_LEVEL OFF)
endif()

option(MUSPINDYNAMIC_BUILD_APPS "Build standalone spin_dynamic_decay CLI application" ${MUSPINDYNAMIC_IS_TOP_LEVEL})
option(MUSPINDYNAMIC_BUILD_TESTING "Build and register unit tests with CTest" ${MUSPINDYNAMIC_IS_TOP_LEVEL})
option(MUSPINDYNAMIC_ENABLE_ROOT_IO "Define MuSpinDynamic::RootIO target when ROOT RIO/Tree/Physics are available" ON)

find_package(ROOT REQUIRED COMPONENTS GenVector Matrix)

# Pure single-event physics core (header-only; no Geant4, Mustard, TFile/TTree, or TRandom3 dependency)
add_library(MuSpinDynamic_Core INTERFACE)
add_library(MuSpinDynamic::Core ALIAS MuSpinDynamic_Core)
add_library(MuSpinDynamic::MuSpinDynamic ALIAS MuSpinDynamic_Core)

target_compile_features(MuSpinDynamic_Core INTERFACE cxx_std_17)
target_include_directories(MuSpinDynamic_Core
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(MuSpinDynamic_Core
    INTERFACE
        ROOT::GenVector
        ROOT::Matrix
)

set(MUSPINDYNAMIC_EXPORT_TARGETS MuSpinDynamic_Core)

# Optional ROOT TTree batch I/O header target (Mode B)
set(MUSPINDYNAMIC_HAS_ROOT_IO OFF)
if(MUSPINDYNAMIC_ENABLE_ROOT_IO OR MUSPINDYNAMIC_BUILD_APPS OR MUSPINDYNAMIC_BUILD_TESTING)
    find_package(ROOT QUIET OPTIONAL_COMPONENTS RIO Tree Physics MathCore)
    if(TARGET ROOT::RIO AND TARGET ROOT::Tree AND TARGET ROOT::Physics)
        add_library(MuSpinDynamic_RootIO INTERFACE)
        add_library(MuSpinDynamic::RootIO ALIAS MuSpinDynamic_RootIO)
        target_link_libraries(MuSpinDynamic_RootIO
            INTERFACE
                MuSpinDynamic_Core
                ROOT::RIO
                ROOT::Tree
                ROOT::Physics
        )
        list(APPEND MUSPINDYNAMIC_EXPORT_TARGETS MuSpinDynamic_RootIO)
        set(MUSPINDYNAMIC_HAS_ROOT_IO ON)
    endif()
endif()

function(muspindynamic_apply_strict_warnings target_name)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${target_name} PRIVATE -Wall -Wextra -Werror)
    endif()
endfunction()

if(MUSPINDYNAMIC_BUILD_APPS)
    if(NOT MUSPINDYNAMIC_HAS_ROOT_IO OR NOT TARGET ROOT::MathCore)
        message(FATAL_ERROR "MUSPINDYNAMIC_BUILD_APPS requires ROOT components RIO, Tree, Physics, and MathCore")
    endif()
    add_executable(spin_dynamic_decay apps/spin_dynamic_decay.cxx)
    target_link_libraries(spin_dynamic_decay PRIVATE MuSpinDynamic::RootIO ROOT::MathCore)
    muspindynamic_apply_strict_warnings(spin_dynamic_decay)
endif()

if(MUSPINDYNAMIC_BUILD_TESTING)
    include(CTest)

    # Core-only isolation test: links ONLY MuSpinDynamic::Core (no ROOT::RIO/Tree/MathCore/TRandom3)
    add_executable(test_core_only tests/test_core_only.cxx)
    target_link_libraries(test_core_only PRIVATE MuSpinDynamic::Core)
    muspindynamic_apply_strict_warnings(test_core_only)
    add_test(NAME test_core_only COMMAND test_core_only)

    # Full end-to-end physics + ROOT I/O test suite
    if(NOT MUSPINDYNAMIC_HAS_ROOT_IO OR NOT TARGET ROOT::MathCore)
        message(FATAL_ERROR "MUSPINDYNAMIC_BUILD_TESTING requires ROOT components RIO, Tree, Physics, and MathCore")
    endif()
    add_executable(test_spin_dynamic_decay tests/test_spin_dynamic_decay.cxx)
    target_link_libraries(test_spin_dynamic_decay PRIVATE MuSpinDynamic::RootIO ROOT::MathCore)
    muspindynamic_apply_strict_warnings(test_spin_dynamic_decay)
    add_test(NAME test_spin_dynamic_decay COMMAND test_spin_dynamic_decay)
endif()

# Install & CMake package export
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS ${MUSPINDYNAMIC_EXPORT_TARGETS} EXPORT MuSpinDynamicTargets)
install(EXPORT MuSpinDynamicTargets
    FILE MuSpinDynamicTargets.cmake
    NAMESPACE MuSpinDynamic::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MuSpinDynamic
)

configure_package_config_file(
    cmake/MuSpinDynamicConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MuSpinDynamicConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MuSpinDynamic
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/MuSpinDynamicConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MuSpinDynamicConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MuSpinDynamicConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MuSpinDynamic
)
