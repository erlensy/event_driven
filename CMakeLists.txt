cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(event_driven)

include_directories(src/*.h)

file(GLOB SOURCES event_driven src/*.cpp src/*.h src/*.py)

add_executable(event_driven ${SOURCES})
