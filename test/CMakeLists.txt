
find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
if(Boost_FOUND)

include_directories(${Boost_INCLUDE_DIRS}/boost)

add_definitions(-DBOOST_TEST_DYN_LINK -g)
add_executable(Test test1.cpp)

target_link_libraries(Test mini_pwa ${Boost_FILESYSTEM_LIBRARY}
                               ${Boost_SYSTEM_LIBRARY}
                               ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

enable_testing()
add_test(NAME mytest COMMAND Test)

endif(Boost_FOUND)
