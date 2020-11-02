cmake_minimum_required(VERSION 3.16)

# General Settings:
set(PROJECT_NAME Encrypted-Estimation)
set(BUILD_SHARED_LIBS YES)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_STANDARD 17)

project(${PROJECT_NAME})
add_definitions("-std=c++17")
enable_testing()

# Set DIR Variables:
set(CRYPTO_DIR ${CMAKE_SOURCE_DIR}/crypto)
set(MATH_DIR ${CMAKE_SOURCE_DIR}/crypto/math)
set(REPRESENTATIONS_DIR ${CMAKE_SOURCE_DIR}/representations)
set(MATLAB_DIR /usr/local/MATLAB/R2020b)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/tests)
set(MATLAB_INCLUDE_DIR ${MATLAB_DIR}/extern/include)
set(MATLAB_LIB_DIR ${MATLAB_DIR}/bin/glnxa64)

# Packages:
find_package(Threads REQUIRED)

# Include Dirs:
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CRYPTO_DIR})
include_directories(${REPRESENTATIONS_DIR})
include_directories(${TEST_DIR})
include_directories(${MATLAB_INCLUDE_DIR})

# Lib Dirs:
link_directories(${MATLAB_LIB_DIR})

# Conan Packages:
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# Sources:
file(GLOB_RECURSE CRYPTO_SRC ${CRYPTO_DIR}/*.cc ${MATH_DIR}/*.cc)
file(GLOB_RECURSE REPRESENTATION_SRC ${REPRESENTATIONS_DIR}/*.cpp)
file(GLOB_RECURSE TEST_SRC ${TEST_DIR}/*Tests.cpp)

message("***[INFO] Sources:")
message(${CRYPTO_SRC})
message(${REPRESENTATION_SRC})
message(${TEST_SRC})

# Targets:
add_executable(ZonoEntities ZonoEntities.cpp ${CRYPTO_SRC} ${REPRESENTATION_SRC})
add_executable(ZonoStrips ZonoStrips.cpp ${CRYPTO_SRC} ${REPRESENTATION_SRC})
add_executable(TestMain ${TEST_DIR}/TestMain.cpp ${CRYPTO_SRC} ${REPRESENTATION_SRC} ${TEST_SRC})

# Linking:
target_link_libraries(ZonoEntities PRIVATE Threads::Threads ${CONAN_LIBS} mat mx)
target_link_libraries(ZonoStrips PRIVATE Threads::Threads ${CONAN_LIBS} mat mx)
target_link_libraries(TestMain PRIVATE Threads::Threads ${CONAN_LIBS} mat mx)

# Installation:
install(TARGETS ZonoEntities ZonoStrips TestMain
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

