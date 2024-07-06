cmake_minimum_required(VERSION 3.15)
project(tsp VERSION 0.0.3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_BUILD_TYPE Release)
if (NOT CMAKE_BUILD_TYPE)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS "-Wall -O3")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-Wall -O3")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W2")
    endif()
endif()

add_library(base_lib OBJECT
    src/util.cpp
    src/eax.cpp
    src/clk.cpp
)

set(binary_name tsp)
add_executable(${binary_name} src/${binary_name}.cpp)
target_link_libraries(${binary_name} base_lib)


set(instances_path ${CMAKE_SOURCE_DIR}/instances)

add_custom_target(test-2opt
    COMMAND ${binary_name} 2opt ${instances_path}/fnl4461.tsp
)

add_custom_target(test-eax
    COMMAND ${binary_name} eax ${instances_path}/eil101.tsp
)

add_custom_target(test-clk
    COMMAND ${binary_name} clk ${instances_path}/pr2392.tsp
)

add_custom_target(test-eax-heavy
    COMMAND ${binary_name} eax ${instances_path}/fnl4461.tsp
)
