add_library(plant
    STATIC
    loader.cpp
    plant.cpp)

target_include_directories(plant
    PRIVATE
    ${CMAKE_SOURCE_DIR})
