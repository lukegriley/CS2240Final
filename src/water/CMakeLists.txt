
add_library(waterdiffusion
    STATIC
    plant.cpp
    loader.cpp
    write.cpp)

target_include_directories(waterdiffusion
    PRIVATE
    ${CMAKE_SOURCE_DIR})

target_link_libraries(waterdiffusion
    PRIVATE
    plant)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(waterdiffusion PUBLIC OpenMP::OpenMP_CXX)
endif()
