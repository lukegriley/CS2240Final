
add_library(stiffrod
    STATIC
    direct.cpp
    rod.cpp
    rodutils.cpp)

target_include_directories(stiffrod
    PRIVATE
    ${CMAKE_SOURCE_DIR})

target_link_libraries(stiffrod
    PRIVATE
    plant)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(stiffrod PUBLIC OpenMP::OpenMP_CXX)
endif()
