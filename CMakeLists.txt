cmake_minimum_required(VERSION 3.16)
project(SPKReader)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Qt6 REQUIRED COMPONENTS Core)
find_library(CFITSIO_LIBRARY cfitsio PATHS /opt/homebrew/lib)
find_library(OPENCV_CORE_LIBRARY opencv_core PATHS /opt/homebrew/lib) 
find_library(OPENCV_AFFINE_LIBRARY opencv_imgproc PATHS /opt/homebrew/lib) 

# Add executable
add_executable(comet_stacker
    comet_stacker.cpp
)

# Link libraries
target_link_libraries(comet_stacker
    Qt6::Core
    ${CFITSIO_LIBRARY}
    ${OPENCV_CORE_LIBRARY}
    ${OPENCV_AFFINE_LIBRARY}
)

target_compile_options(comet_stacker PRIVATE
    -Wno-dangling-else
    -Wno-extern-c-compat
    -Wno-asm-operand-widths
    -O3 -g
)

# Platform-specific settings
if(APPLE)
    # macOS might need additional flags
    target_include_directories(comet_stacker PUBLIC
        /opt/homebrew/Cellar/opencv/4.12.0_4/include/opencv4)
    target_link_libraries(comet_stacker "-framework CoreFoundation")
endif()
