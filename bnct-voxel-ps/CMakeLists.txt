cmake_minimum_required(VERSION 3.16...3.27)
project(BNCTVoxelPhaseSpace)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and visualization drivers" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(bnctVoxelPS
  bnctVoxelPS.cc
  ${sources}
  ${headers}
)

target_include_directories(bnctVoxelPS PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)
if(WITH_GEANT4_UIVIS)
  target_compile_definitions(bnctVoxelPS PRIVATE BNCT_WITH_GEANT4_UIVIS=1)
endif()
target_link_libraries(bnctVoxelPS ${Geant4_LIBRARIES})
target_compile_features(bnctVoxelPS PRIVATE cxx_std_17)

file(GLOB macro_files ${PROJECT_SOURCE_DIR}/macros/*.mac)
foreach(_script ${macro_files})
  configure_file(${_script} ${PROJECT_BINARY_DIR}/. COPYONLY)
endforeach()

install(TARGETS bnctVoxelPS DESTINATION bin)
