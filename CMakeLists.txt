cmake_minimum_required(VERSION 3.16)

project(comet_stacker VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -g") 
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -g")

# Find Qt6
find_package(Qt6 REQUIRED COMPONENTS Core Widgets Network Charts)

# Set up Qt6 automoc
set(CMAKE_AUTOMOC ON)

# PCL Configuration
set(PCL_ROOT "/Users/jonathan/PCL" CACHE PATH "PCL root directory")
set(PCL_INCLUDE_DIRS 
    "${PCL_ROOT}/include" 
    "${PCL_ROOT}/src/3rdparty"
    "${PCL_ROOT}/src/modules/file-formats"
    "${PCL_ROOT}/src/modules/file-formats/FITS"
    "${PCL_ROOT}/src"
)
set(PCL_LIBRARY_DIR "${PCL_ROOT}/src/pcl/macosx/g++/arm64/Release")

# Find PCL library
find_library(PCL_LIBRARY 
    NAMES PCL-pxi libPCL-pxi.a
    PATHS ${PCL_LIBRARY_DIR}
    NO_DEFAULT_PATH
)

if(NOT PCL_LIBRARY)
    message(FATAL_ERROR "PCL library not found in ${PCL_LIBRARY_DIR}")
endif()

message(STATUS "Found PCL library: ${PCL_LIBRARY}")

# Create SHA library explicitly
message(STATUS "Setting up PCL SHA sources...")

set(SHA1_FILE "${PCL_ROOT}/src/3rdparty/RFC6234/sha1.c")
set(SHA224_256_FILE "${PCL_ROOT}/src/3rdparty/RFC6234/sha224-256.c") 
set(SHA384_512_FILE "${PCL_ROOT}/src/3rdparty/RFC6234/sha384-512.c")

foreach(sha_file ${SHA1_FILE} ${SHA224_256_FILE} ${SHA384_512_FILE})
    if(NOT EXISTS ${sha_file})
        message(FATAL_ERROR "SHA source not found: ${sha_file}")
    else()
        message(STATUS "Found SHA source: ${sha_file}")
    endif()
endforeach()

# Create static library for SHA functions
add_library(pcl_sha STATIC
    ${SHA1_FILE}
    ${SHA224_256_FILE}
    ${SHA384_512_FILE}
)

set_target_properties(pcl_sha PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(pcl_sha PUBLIC
    "${PCL_ROOT}/src/3rdparty/RFC6234"
    "${PCL_ROOT}/src/3rdparty"
)

# Create CMINPACK library from PCL sources
message(STATUS "Setting up CMINPACK sources...")

set(CMINPACK_DIR "${PCL_ROOT}/src/3rdparty/cminpack")
set(CMINPACK_SOURCES
    "${CMINPACK_DIR}/lmdif1.c"
    "${CMINPACK_DIR}/lmdif.c"
    "${CMINPACK_DIR}/lmpar.c"
    "${CMINPACK_DIR}/qrfac.c"
    "${CMINPACK_DIR}/qrsolv.c"
    "${CMINPACK_DIR}/enorm.c"
    "${CMINPACK_DIR}/dpmpar.c"
    "${CMINPACK_DIR}/lmdif1_.c"
    "${CMINPACK_DIR}/lmdif_.c"
    "${CMINPACK_DIR}/lmpar_.c"
    "${CMINPACK_DIR}/qrfac_.c"
    "${CMINPACK_DIR}/qrsolv_.c"
    "${CMINPACK_DIR}/enorm_.c"
    "${CMINPACK_DIR}/dpmpar_.c"
    "${CMINPACK_DIR}/fdjac1.c"
    "${CMINPACK_DIR}/fdjac1_.c"
    "${CMINPACK_DIR}/fdjac2.c"
    "${CMINPACK_DIR}/fdjac2_.c"
)

# Check which CMINPACK sources exist
set(AVAILABLE_CMINPACK_SOURCES)
foreach(src_file ${CMINPACK_SOURCES})
    if(EXISTS ${src_file})
        list(APPEND AVAILABLE_CMINPACK_SOURCES ${src_file})
        message(STATUS "Found CMINPACK source: ${src_file}")
    else()
        message(STATUS "CMINPACK source not found: ${src_file}")
    endif()
endforeach()

if(AVAILABLE_CMINPACK_SOURCES)
    message(STATUS "Creating CMINPACK library with ${CMAKE_CURRENT_LIST_LENGTH} sources")
    
    add_library(pcl_cminpack STATIC ${AVAILABLE_CMINPACK_SOURCES})
    
    set_target_properties(pcl_cminpack PROPERTIES
        C_STANDARD 99
        C_STANDARD_REQUIRED ON
        POSITION_INDEPENDENT_CODE ON
    )
    
    target_include_directories(pcl_cminpack PUBLIC
        "${CMINPACK_DIR}"
        "${PCL_ROOT}/src/3rdparty"
    )
    
    # Compiler flags for CMINPACK
    target_compile_options(pcl_cminpack PRIVATE
        -O3 -g
        -Wno-unused-parameter
        -Wno-unused-variable
        -Wno-unused-but-set-variable
        -Wno-uninitialized
    )
    
    # Define preprocessor macros for CMINPACK
    target_compile_definitions(pcl_cminpack PRIVATE
        __CMINPACK_NO_DLL__
        CMINPACK_NO_DLL
    )
    
    set(HAVE_CMINPACK TRUE)
else()
    message(WARNING "No CMINPACK sources found in ${CMINPACK_DIR}")
    set(HAVE_CMINPACK FALSE)
endif()

# Add PCL Mock Library
add_library(pcl_mock STATIC
    PCLMockAPI.cpp
    PCLThreadMock.cpp
)

set_target_properties(pcl_mock PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(pcl_mock PUBLIC
    ${PCL_INCLUDE_DIRS}
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_compile_definitions(pcl_mock PRIVATE
    __PCL_MACOSX
    __PCL_NO_PERFORMANCE_CRITICAL_MATH_ROUTINES
)

# Compiler flags
target_compile_options(pcl_mock PRIVATE
    -Wno-dangling-else
    -Wno-extern-c-compat
    -Wno-asm-operand-widths
    -O3 -g
)

# Third-party libraries
find_library(LZ4_LIBRARY lz4 PATHS /opt/homebrew/lib)
find_library(ZSTD_LIBRARY zstd PATHS /opt/homebrew/lib)
find_library(LCMS2_LIBRARY lcms2 PATHS /opt/homebrew/lib)
find_library(CRYPTO_LIBRARY crypto PATHS /opt/homebrew/lib)
find_library(CFITSIO_LIBRARY cfitsio PATHS /opt/homebrew/lib)
find_library(TIFF_LIBRARY NAMES tiff PATHS /opt/homebrew/lib)
find_library(CURL_LIBRARY NAMES curl PATHS /opt/homebrew/lib)
find_library(OPENCV_CORE_LIBRARY opencv_core PATHS /opt/homebrew/lib) 
find_library(OPENCV_AFFINE_LIBRARY opencv_imgproc PATHS /opt/homebrew/lib) 

# Check for required libraries
foreach(lib LZ4_LIBRARY ZSTD_LIBRARY LCMS2_LIBRARY CRYPTO_LIBRARY TIFF_LIBRARY CURL_LIBRARY)
    if(NOT ${lib})
        message(FATAL_ERROR "${lib} not found. Install with: brew install lz4 zstd little-cms2 openssl libtiff curl")
    endif()
endforeach()

# CFITSIO is optional but recommended for FITS support
if(NOT CFITSIO_LIBRARY)
    message(WARNING "CFITSIO library not found. Install with: brew install cfitsio")
endif()

# Include directories
set(PROJECT_INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

set(SOURCES
    comet_stacker.cpp
    GaiaGDR3Catalog.cpp
    StarCatalogValidator.cpp
    FITS.cpp
    ${PCL_ROOT}/src/modules/file-formats/TIFF/TIFF.cpp
    ${PCL_ROOT}/src/modules/file-formats/TIFF/TIFFFormat.cpp
    ${PCL_ROOT}/src/modules/file-formats/TIFF/TIFFInstance.cpp
    ${PCL_ROOT}/src/modules/file-formats/TIFF/TIFFModule.cpp
    ${PCL_ROOT}/src/modules/file-formats/TIFF/TIFFOptionsDialog.cpp
    ${PCL_ROOT}/src/modules/file-formats/TIFF/TIFFPreferencesDialog.cpp
    ${PCL_ROOT}/src/modules/file-formats/TIFF/TIFFRangeOptionsDialog.cpp
)

set(HEADERS
    GaiaGDR3Catalog.h
    StarCatalogValidator.h
)

add_executable(comet_stacker ${SOURCES} ${HEADERS})

# Set target properties
set_target_properties(comet_stacker PROPERTIES
    MACOSX_BUNDLE TRUE
    MACOSX_BUNDLE_GUI_IDENTIFIER com.stellina.xisfcreator
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION}
    WIN32_EXECUTABLE TRUE
)

# Include directories
target_include_directories(comet_stacker PRIVATE
    ${PCL_INCLUDE_DIRS}
    ${PROJECT_INCLUDE_DIRS}
    /opt/homebrew/include
    /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3
    /opt/homebrew/Cellar/opencv/4.12.0_4/include/opencv4
    /opt/homebrew/Cellar/qt/6.9.0/lib/QtNetwork.framework/Versions/A/Headers
)

set(ASTROMETRY_PREFIX /opt/homebrew) 
set(HOMEBREW_PREFIX /opt/homebrew) 

# Link libraries
target_link_libraries(comet_stacker PRIVATE
    Qt6::Core
    Qt6::Widgets
    Qt6::Network
    Qt6::Charts
    pcl_mock          # Add mock library
    ${PCL_LIBRARY}
    pcl_sha
    ${LZ4_LIBRARY}
    ${ZSTD_LIBRARY} 
    ${LCMS2_LIBRARY}
    ${CRYPTO_LIBRARY}
    ${TIFF_LIBRARY}
    ${CURL_LIBRARY}
    ${CFITSIO_LIBRARY}
    ${OPENCV_CORE_LIBRARY}
    ${OPENCV_AFFINE_LIBRARY}
    ${ASTROMETRY_PREFIX}/lib/libastrometry.a
    ${ASTROMETRY_PREFIX}/lib/libanfiles.a
    ${ASTROMETRY_PREFIX}/lib/libkd.a
    ${ASTROMETRY_PREFIX}/lib/libanutils.a
    ${ASTROMETRY_PREFIX}/lib/libanbase.a
    ${ASTROMETRY_PREFIX}/lib/libqfits.a
    ${HOMEBREW_PREFIX}/lib/libgsl.a
    ${HOMEBREW_PREFIX}/lib/libgslcblas.a
    ${HOMEBREW_PREFIX}/lib/libwcs.a
    
    z
)

# Link CMINPACK if available
if(HAVE_CMINPACK)
    target_link_libraries(comet_stacker PRIVATE pcl_cminpack)
    target_compile_definitions(comet_stacker PRIVATE HAVE_CMINPACK)
    message(STATUS "CMINPACK support enabled - PSF fitting available")
else()
    target_compile_definitions(comet_stacker PRIVATE NO_CMINPACK)
    message(STATUS "CMINPACK not available - PSF fitting disabled")
endif()

# Add CFITSIO if available
if(CFITSIO_LIBRARY)
    target_link_libraries(comet_stacker PRIVATE ${CFITSIO_LIBRARY})
endif()

# macOS-specific frameworks (required for FITS support)
if(APPLE)
    target_link_libraries(comet_stacker PRIVATE
        "-framework CoreFoundation"
        "-framework Foundation"
    )
endif()

# Compiler definitions
target_compile_definitions(comet_stacker PRIVATE
    __PCL_MACOSX
    __PCL_NO_PERFORMANCE_CRITICAL_MATH_ROUTINES
)

# Compiler flags
target_compile_options(comet_stacker PRIVATE
    -Wno-dangling-else
    -Wno-extern-c-compat
    -Wno-asm-operand-widths
    -O0 -g
)

# Debug output
message(STATUS "PCL Include Dirs: ${PCL_INCLUDE_DIRS}")
message(STATUS "PCL Library: ${PCL_LIBRARY}")
message(STATUS "Qt6 Version: ${Qt6_VERSION}")
message(STATUS "CMINPACK Available: ${HAVE_CMINPACK}")
