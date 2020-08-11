# 'Google_test' is the subproject name
project(Google_tests)

# 'lib' is the folder with Google Test sources
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories("${PROJECT_SOURCE_DIR}/src/lib")
find_package(OpenSSL)
# 'Google_Tests_run' is the target name
# 'test1.cpp tests2.cpp' are source files with tests
add_executable(Google_Tests_run ../tests/gf16-test.cc
        ../src/lib/gf16.c
        ../src/lib/verify.c
        ../src/main/timing-functions.c
        ../src/lib/keygen.c
        ../src/lib/rng.c
        ../src/lib/utils_hash.c
        ../src/lib/utils_prng.c
        ../tests/keygen_test.cpp ../src/lib/sign.c ../src/lib/sign.h ../tests/sign_test.cpp ../tests/verify_test.cpp)
target_link_libraries(Google_Tests_run gtest gtest_main OpenSSL::Crypto)