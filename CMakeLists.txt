# CMake 最低版本号要求

cmake_minimum_required (VERSION 3.10)

# 项目信息

project(mini_pwa)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# add_compile_options(-fpermissive)
add_subdirectory(3rd/Minuit2)
add_subdirectory(3rd/cnpy)

include_directories(include)



# 指定生成目标
add_library(mini_pwa src/cg_coeff.cpp src/decay.cpp src/dfunction.cpp src/variable.cpp src/fcn.cpp)
add_executable(random_fit src/main.cpp)

target_include_directories(random_fit PUBLIC 3rd/cnpy)
target_link_libraries(mini_pwa PUBLIC Minuit2::Minuit2)
target_link_libraries(random_fit PUBLIC mini_pwa cnpy)

add_subdirectory(test)
