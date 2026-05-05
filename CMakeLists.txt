cmake_minimum_required(VERSION 3.16)

project(X17_ATOMKI_simulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

add_executable(x17sim
    x17sim.cc
    src/DetectorConstruction.cc
    src/ActionInitialization.cc
    src/PrimaryGeneratorAction.cc
    src/RunAction.cc
    src/EventAction.cc
    src/SteppingAction.cc
)

target_include_directories(x17sim PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(x17sim PRIVATE
    ${Geant4_LIBRARIES}
)

install(TARGETS x17sim DESTINATION bin)
