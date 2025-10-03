cmake_minimum_required(VERSION 3.15)
project(ZPinchStabilitySim VERSION 1.0.0 LANGUAGES CXX)

# ==============================================================================
# COMPILER SETTINGS - MSVC COMPATIBLE
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Global compiler flags (applied to all targets)
if(MSVC)
    # MSVC-specific settings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /MP /EHsc")
    set(CMAKE_CXX_FLAGS_DEBUG "/DDEBUG /Zi /Od /RTC1")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/std:c++17>")
else()
    # GCC/Clang settings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -O2")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

# ==============================================================================
# PROJECT STRUCTURE
# ==============================================================================
set(INCLUDE_DIR include)
set(SOURCE_DIR src)
set(PYTHON_SCRIPTS_DIR pyscripts)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================

# Find required packages
find_package(Threads REQUIRED)

# Check for optional dependencies
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: enabled")
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
else()
    message(WARNING "OpenMP not found: parallelization disabled")
endif()

# ==============================================================================
# SOURCE FILES
# ==============================================================================

# Main executable sources
set(MAIN_SOURCES
    ${SOURCE_DIR}/main.cpp
    ${SOURCE_DIR}/equilibrium_solver.cpp
    ${SOURCE_DIR}/linear_stability.cpp  
    ${SOURCE_DIR}/nonlinear_evolution.cpp
    ${SOURCE_DIR}/diagnostics.cpp
)

# Headers
set(HEADER_FILES
    ${INCLUDE_DIR}/zpinch_params.hh
    ${INCLUDE_DIR}/equilibrium_solver.hh
    ${INCLUDE_DIR}/linear_stability.hh
    ${INCLUDE_DIR}/nonlinear_evolution.hh
    ${INCLUDE_DIR}/diagnostics.hh
)

# ==============================================================================
# EXECUTABLE TARGET
# ==============================================================================

# Create the main executable
add_executable(zpinch_sim ${MAIN_SOURCES} ${HEADER_FILES})

# Include directories
target_include_directories(zpinch_sim PRIVATE ${INCLUDE_DIR})

# Link libraries
target_link_libraries(zpinch_sim PRIVATE Threads::Threads)

if(OpenMP_CXX_FOUND)
    target_link_libraries(zpinch_sim PRIVATE OpenMP::OpenMP_CXX)
endif()

# Math library (may be needed on some systems)
if(UNIX AND NOT APPLE)
    target_link_libraries(zpinch_sim PRIVATE m)
endif()

# ==============================================================================
# TARGET-SPECIFIC COMPILER OPTIONS (AFTER TARGET CREATION)
# ==============================================================================
if(MSVC)
    # MSVC-specific warnings suppression
    target_compile_options(zpinch_sim PRIVATE 
        /wd4996  # Disable deprecated function warnings
        /wd4267  # Disable size_t to int conversion warnings
        /wd4244  # Disable possible loss of data warnings
        /wd4189  # Disable unused variable warnings
        /wd4100  # Disable unused parameter warnings
    )
    
    # MSVC-specific definitions
    target_compile_definitions(zpinch_sim PRIVATE 
        _USE_MATH_DEFINES
        NOMINMAX  # Disable min/max macros
    )
endif()

# ==============================================================================
# INSTALLATION
# ==============================================================================

# Install executable
install(TARGETS zpinch_sim
    RUNTIME DESTINATION bin
)

# Install headers
install(DIRECTORY ${INCLUDE_DIR}/
    DESTINATION include/zpinch
    FILES_MATCHING PATTERN "*.hh"
)

# ==============================================================================
# TESTING (SIMPLIFIED FOR MSVC)
# ==============================================================================

# Enable testing if requested
option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    
    # Simple test to check if executable runs
    add_test(NAME BasicExecution
        COMMAND zpinch_sim help
    )
endif()

add_custom_command(TARGET zpinch_sim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:zpinch_sim>/data
)
add_custom_command(TARGET zpinch_sim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:zpinch_sim>/data/snap
)
# ==============================================================================
# CUSTOM TARGETS
# ==============================================================================

# Target to run the full analysis pipeline
add_custom_target(run_full
    COMMAND zpinch_sim full
    DEPENDS zpinch_sim
    COMMENT "Running full Z-pinch analysis pipeline"
)

# ==============================================================================
# POST-BUILD CHECKS
# ==============================================================================

# Check if the binary was built successfully
add_custom_command(TARGET zpinch_sim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "ZPinchStabilitySim built successfully!"
    COMMAND ${CMAKE_COMMAND} -E echo "Usage: zpinch_sim [equilibrium|linear|nonlinear|full|help]"
)

# ==============================================================================
# OUTPUT MESSAGES
# ==============================================================================

message(STATUS "")
message(STATUS "===================================================")
message(STATUS "ZPinchStabilitySim Configuration Summary")
message(STATUS "===================================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ standard: 17")
message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "Build tests: ${BUILD_TESTS}")
message(STATUS "===================================================")
message(STATUS "")
message(STATUS "Build the project with: cmake --build .")
message(STATUS "Run full analysis: cmake --build . --target run_full")
message(STATUS "")
