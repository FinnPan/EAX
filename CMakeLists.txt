cmake_minimum_required(VERSION 3.15)
project(tsp VERSION 0.0.2 LANGUAGES CXX)


if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(my_cxxflags "-Wall")
    set(extra_libs m)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(my_cxxflags "/W2")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${my_cxxflags}")


add_library(base_lib OBJECT
    src/util.cpp
    src/eax.cpp
    src/clk.cpp
)
link_libraries(${extra_libs} base_lib)


set(binary_name tsp)
add_executable(${binary_name} src/${binary_name}.cpp)


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
