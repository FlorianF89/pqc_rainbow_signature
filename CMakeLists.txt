cmake_minimum_required(VERSION 3.15)
project(pqc_rainbow C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 14)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-tree-vectorize -Ofast")
    message("c flags = ${CMAKE_C_FLAGS}")
endif ()

find_package(OpenSSL REQUIRED)
set(OPENSSL_INCLUDE_DIR "/usr/local/include/")
include_directories(${OPENSSL_INCLUDE_DIR})

add_library(rainbow src/lib/gf16.c src/lib/gf16.h src/lib/parameters.h src/lib/keygen.c src/lib/keygen.h
            src/lib/sign.h src/lib/sign.c src/lib/rng.c src/lib/rng.h src/lib/utils_hash.h src/lib/utils_hash.c
            src/lib/verify.h src/lib/verify.c)
add_executable(test_pqc_rainbow src/main/main_test.c src/main/test_functions.c)

add_executable(keygen_pqc_rainbow src/main/main_keygen.c)

add_dependencies(rainbow OpenSSL::SSL OpenSSL::Crypto)
target_link_libraries(rainbow OpenSSL::SSL OpenSSL::Crypto)
add_dependencies(test_pqc_rainbow rainbow)
target_link_libraries(test_pqc_rainbow rainbow)

add_dependencies(keygen_pqc_rainbow rainbow)
target_link_libraries(keygen_pqc_rainbow rainbow)