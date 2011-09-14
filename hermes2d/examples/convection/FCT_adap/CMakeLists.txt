project(FCT_adap)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp lumped_projection.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

