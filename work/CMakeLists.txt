project(work)

set(includes ./ ../non_linear_solvers/)

add_executable(work main.cpp)
target_include_directories(work PUBLIC ${includes})
target_link_libraries(work non_linear_solvers)
