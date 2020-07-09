cmake_minimum_required(VERSION 3.15)
project(pqc_rainbow C)

set(CMAKE_C_STANDARD 99)

add_library(rainbow src/lib/gf16.c src/lib/gf16.h)
add_executable(pqc_rainbow src/main/timing-functions.c src/main/timing-functions.h src/main/signature_system.c)
add_dependencies(pqc_rainbow rainbow)
target_link_libraries(pqc_rainbow rainbow)