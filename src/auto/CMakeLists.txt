# Sets C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specifies .cpp and .h files to be passed to the compiler
add_executable(run-simulation
    config.cpp
    export.cpp
    main.cpp
    simulation.cpp
)

target_include_directories(run-simulation
    PRIVATE
    ${CMAKE_SOURCE_DIR})

# Specifies required Qt components
find_package(Qt6 REQUIRED COMPONENTS Core)
find_package(Qt6 REQUIRED COMPONENTS Concurrent)
find_package(Qt6 REQUIRED COMPONENTS Xml)

# Specifies libraries to be linked (Qt components, glew, etc)
target_link_libraries(run-simulation PRIVATE
    Qt::Concurrent
    Qt::Core
    Qt::Xml
    waterdiffusion
    stiffrod
)

