cmake_minimum_required(VERSION 3.12.0)

project(epsic, LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)


file(GLOB_RECURSE source CONFIGURE_DEPENDS include/*.h src/*.cpp src/*.c)
add_library(epsiclib SHARED ${source})

target_include_directories(epsiclib PUBLIC include)

add_executable(epsic epsic.cpp)
target_link_libraries(epsic PUBLIC epsiclib)

add_subdirectory(test)
set_target_properties(epsiclib PROPERTIES OUTPUT_NAME "epsic")

