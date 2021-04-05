# CMake 最低版本号要求

cmake_minimum_required (VERSION 2.8)

# 项目信息

project(mini_pwa)
set(CMAKE_BUILD_TYPE "Debug")

add_compile_options(-fpermissive)

include_directories(include)

add_subdirectory(3rd/Minuit2)

# 指定生成目标

add_executable(random_fit src/main.cpp src/cg_coeff.cpp src/decay.cpp src/dfunction.cpp src/variable.cpp)
target_link_libraries(random_fit PUBLIC Minuit2::Minuit2)
