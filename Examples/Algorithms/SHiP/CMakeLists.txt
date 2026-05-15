# Use the ACTS-specific library macro
acts_add_library(
    ExamplesSHiP
    SHARED
    src/SeedToParams.cpp
    src/SHiPMeasurementProvider.cpp
    src/SHiPSeeder.cpp
    src/RecoTrack.cpp
    src/RecoVertex.cpp
)

# This links the required ACTS components
target_link_libraries(
    ActsExamplesSHiP
    PUBLIC
        Acts::Core
        ActsExamplesFramework
        ActsExamplesDigitization
    PRIVATE 
        ROOT::Core
)

# Crucial: Ensure the include directory is visible to other targets
target_include_directories(
    ActsExamplesSHiP 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

root_generate_dictionary(
    ActsExamplesSHiPDict 
    include/ActsExamples/SHiP/RecoTrack.hpp
    include/ActsExamples/SHiP/RecoVertex.hpp
    MODULE ActsExamplesSHiP 
    LINKDEF LinkDef.h
)

acts_compile_headers(ExamplesSHiP GLOB "include/**/*.hpp")
# Add standard Acts cleanup properties
set_target_properties(ActsExamplesSHiPDict PROPERTIES CXX_CLANG_TIDY "")
target_compile_options(ActsExamplesSHiPDict PRIVATE "-Wno-old-style-cast")
