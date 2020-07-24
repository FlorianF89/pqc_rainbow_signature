cmake_minimum_required(VERSION 3.15)
project(pqc_rainbow C)

set(CMAKE_C_STANDARD 99)

if (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -o3")
    message("c flags = ${CMAKE_C_FLAGS}")
endif ()

find_package(OpenSSL REQUIRED)
set(OPENSSL_INCLUDE_DIR "/usr/local/include/")
include_directories(${OPENSSL_INCLUDE_DIR})

add_library(rainbow src/lib/gf16.c src/lib/gf16.h src/lib/parameters.h src/lib/keygen.c src/lib/keygen.h src/lib/rng.c src/lib/rng.h src/lib/utils_prng.c src/lib/utils_prng.h src/lib/utils_hash.h src/lib/utils_hash.c src/lib/hash_len_config.h src/lib/error_codes.h)
add_executable(pqc_rainbow src/main/timing-functions.c src/main/timing-functions.h src/main/signature_system.c)
add_dependencies(pqc_rainbow rainbow)
target_link_libraries(pqc_rainbow rainbow)
add_subdirectory(Google_tests)