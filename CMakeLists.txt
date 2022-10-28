# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 3.20.1)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

if(EXISTS "/Users/hauke/src/itk/InsightToolkit-5.2.1/bin")
  set(ITK_DIR /Users/hauke/src/itk/InsightToolkit-5.2.1/bin)
else()
  set(ITK_DIR /opt/InsightToolkit-5.2.1/bin)
endif()

IF(UNIX AND NOT APPLE)
  set(LINUX TRUE)
  message(STATUS ">>> Linux")
ELSE()
  message(STATUS ">>> MacOS")
ENDIF()

# SET(CMAKE_CXX_FLAGS "-std=c++0x")
IF(LINUX)
  SET(CMAKE_CXX_FLAGS "-std=c++11")
ELSE()
  # find the correct path for isysroot with `xcrun --show-sdk-path`
  SET(CMAKE_CXX_FLAGS "-std=c++11 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -I/usr/local/include")
ENDIF()

project(PR2MASK)

FIND_PACKAGE(ITK)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build ITKApps without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Boost
set(boost_min_ver 1.50.0)
set(boost_libs system filesystem timer date_time)
find_package(Boost ${boost_min_ver})

if(Boost_FOUND)
  find_package(Boost ${boost_min_ver} COMPONENTS ${boost_libs})
endif()

find_package(TBB REQUIRED)

if(TBB_FOUND)
  MESSAGE(STATUS ">>> Found TBB")
ENDIF()

# target_link_libraries(your_target ${Boost_LIBRARIES})
add_executable(pr2mask pr2mask.cxx)
target_link_libraries(pr2mask ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb)
target_compile_features(pr2mask PRIVATE cxx_std_11)

add_executable(mask2rtstruct mask2rtstruct.cxx)
target_link_libraries(mask2rtstruct ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb)
target_compile_features(mask2rtstruct PRIVATE cxx_std_11)