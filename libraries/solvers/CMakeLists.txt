project(solvers)

set(includes includes/ ../utils/)
set(sources # src/common.cpp
            # src/non_linear_solvers.cpp
            src/BisectionMethod.cpp
            src/FalsePositionMethod.cpp
            src/DekkersMethod.cpp
            src/BrentsMethod.cpp
            src/NewtonRaphsonMethod.cpp
            # src/Function.cpp
)

add_library(${PROJECT_NAME} SHARED ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes} ${common_includes_dir})
target_link_libraries(${PROJECT_NAME} PUBLIC utils)

# install(TARGETS operations DESTINATION ${CUSTOM_LIBRARY_DIR_PYTHON})