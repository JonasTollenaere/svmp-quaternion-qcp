project(svmp-quaternion-qcp)

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

### CPP ###
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

### CGAL ###
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE "TRUE")
find_package(CGAL REQUIRED)
link_libraries(CGAL::CGAL)
include_directories(${CGAL_INCLUDE_DIRS})

### MeshCore ###
add_compile_definitions(GLM_ENABLE_EXPERIMENTAL)
set(MESHCORE_BUILD_DEMOS OFF)
include(FetchContent)
FetchContent_Declare(MeshCore
        GIT_REPOSITORY https://github.com/JonasTollenaere/MeshCore.git
        GIT_TAG        v0.1.3
)
FetchContent_MakeAvailable(MeshCore)

### Project headers and source files ###
include_directories(include)
file(GLOB SOURCES src/*.cpp)


#find_package(GUROBI REQUIRED) # FindGUROBI.cmake not found when using vcpkg toolchain, so we pass paths as explicit CMAKE variables
# Examples:
# set(GUROBI_INCLUDE_DIRS ${GUROBI_DIR}/include)
# set(GUROBI_CXX_LIBRARY ${GUROBI_DIR}/lib/gurobi_c++md2017.lib)
# set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_DIR}/lib/gurobi_c++mdd2017.lib)
# set(GUROBI_LIBRARY ${GUROBI_DIR}/lib/gurobi110.lib)
if (DEFINED GUROBI_INCLUDE_DIRS AND DEFINED GUROBI_CXX_LIBRARY AND DEFINED GUROBI_CXX_DEBUG_LIBRARY AND DEFINED GUROBI_LIBRARY)

    ### ConvexQCP target ###
    add_executable(ConvexQuaternionQCP targets/ConvexQCP.cpp ${SOURCES})
    target_include_directories(ConvexQuaternionQCP PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(ConvexQuaternionQCP PRIVATE MeshCore)
    target_link_libraries(ConvexQuaternionQCP PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(ConvexQuaternionQCP PRIVATE ${GUROBI_LIBRARY})

    ### ConvexQCPBenchmark target ###
    add_executable(ConvexQuaternionQCPBenchmark targets/ConvexQCPBenchmark.cpp ${SOURCES})
    target_include_directories(ConvexQuaternionQCPBenchmark PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(ConvexQuaternionQCPBenchmark PRIVATE MeshCore)
    target_link_libraries(ConvexQuaternionQCPBenchmark PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(ConvexQuaternionQCPBenchmark PRIVATE ${GUROBI_LIBRARY})

    ### ConvexQCPInverse target ###
    add_executable(ConvexQuaternionQCPInverse targets/ConvexQCPInverse.cpp ${SOURCES})
    target_include_directories(ConvexQuaternionQCPInverse PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(ConvexQuaternionQCPInverse PRIVATE MeshCore)
    target_link_libraries(ConvexQuaternionQCPInverse PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(ConvexQuaternionQCPInverse PRIVATE ${GUROBI_LIBRARY})

    ### ConvexQCPInverseBenchmark target ###
    add_executable(ConvexQuaternionQCPInverseBenchmark targets/ConvexQCPInverseBenchmark.cpp ${SOURCES})
    target_include_directories(ConvexQuaternionQCPInverseBenchmark PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(ConvexQuaternionQCPInverseBenchmark PRIVATE MeshCore)
    target_link_libraries(ConvexQuaternionQCPInverseBenchmark PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(ConvexQuaternionQCPInverseBenchmark PRIVATE ${GUROBI_LIBRARY})

    ### ConcaveQCP target ###
    add_executable(ConcaveQuaternionQCP targets/ConcaveQCP.cpp ${SOURCES})
    target_include_directories(ConcaveQuaternionQCP PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(ConcaveQuaternionQCP PRIVATE MeshCore)
    target_link_libraries(ConcaveQuaternionQCP PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(ConcaveQuaternionQCP PRIVATE ${GUROBI_LIBRARY})

    ### ConcaveQCPBenchmark target ###
    add_executable(ConcaveQuaternionQCPBenchmark targets/ConcaveQCPBenchmark.cpp ${SOURCES})
    target_include_directories(ConcaveQuaternionQCPBenchmark PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(ConcaveQuaternionQCPBenchmark PRIVATE MeshCore)
    target_link_libraries(ConcaveQuaternionQCPBenchmark PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(ConcaveQuaternionQCPBenchmark PRIVATE ${GUROBI_LIBRARY})

else ()
    message("GUROBI not found. Skipping Strip_Packing_MILP.exe")
endif ()


### Visualize_Symmetry target ###
add_executable(VisualizeSymmetry targets/VisualizeSymmetry.cpp ${SOURCES})
target_link_libraries(VisualizeSymmetry PRIVATE MeshCore)

### Visualize_Concavities target ###
add_executable(VisualizeConcavities targets/VisualizeConcavities.cpp ${SOURCES})
target_link_libraries(VisualizeConcavities PRIVATE MeshCore)