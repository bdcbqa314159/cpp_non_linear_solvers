project(work)

set(includes ./ ../libraries/utils/ ../libraries/solvers/)

add_executable(work main.cpp)
target_include_directories(work PUBLIC ${includes})
target_link_libraries(work utils solvers)

add_executable(my_test test.cpp)

