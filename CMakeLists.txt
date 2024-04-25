cmake_minimum_required(VERSION 3.15)
project(tsp VERSION 0.0.1 LANGUAGES CXX)


if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(FLAGS_COMMON "-Wall -Werror -Wextra")
    set(CMAKE_CXX_FLAGS "${FLAGS_COMMON} -O3 -g0")
    set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")
    set(extra_libs m)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

else ()
    message(FATAL_ERROR "unsupported compiler: " ${CMAKE_CXX_COMPILER_ID})
endif()


add_library(base_lib OBJECT
    src/util.cpp
    src/eax.cpp
    src/clk.cpp
)
link_libraries(${extra_libs} base_lib)


set(target_2opt test_2opt)
add_executable(${target_2opt} src/${target_2opt}.cpp)

set(target_eax test_eax)
add_executable(${target_eax} src/${target_eax}.cpp)

set(target_clk test_clk)
add_executable(${target_clk} src/${target_clk}.cpp)


set(instances_path ${CMAKE_SOURCE_DIR}/instances)

add_custom_target(test-2opt
    COMMAND ${target_2opt} ${instances_path}/fnl4461.tsp
)

add_custom_target(test-eax
    COMMAND ${target_eax} ${instances_path}/eil101.tsp
)

add_custom_target(test-clk
    COMMAND ${target_clk} ${instances_path}/pr2392.tsp
)

add_custom_target(test-all
    COMMAND ${target_2opt} ${instances_path}/fnl4461.tsp
    COMMAND ${target_eax} ${instances_path}/eil101.tsp
    COMMAND ${target_clk} ${instances_path}/pr2392.tsp
)

add_custom_target(test-eax-heavy
    COMMAND ${target_eax} ${instances_path}/fnl4461.tsp
)
