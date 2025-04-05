project(non_linear_solvers)

set(includes includes/)
set(sources src/non_linear_solvers.cpp
            src/BisectionMethod.cpp
            src/FalsePositionMethod.cpp
            src/DekkersMethod.cpp
            src/BrentsMethod.cpp
            src/NewtonRaphsonMethod.cpp
            src/Function.cpp
)

add_library(${PROJECT_NAME} ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes})
target_link_libraries(${PROJECT_NAME} PUBLIC GSL::gsl GSL::gslcblas Boost::boost Boost::date_time)
