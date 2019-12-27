cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0015 NEW)

# set the project name
project(raytracer)

# add the executable
add_executable(raytracer main.cpp)

target_include_directories(raytracer PUBLIC include)

target_link_directories(raytracer PUBLIC lib)
target_link_libraries(raytracer PUBLIC SDL2-2.0.0)
