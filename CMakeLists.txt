cmake_minimum_required(VERSION 3.20.1)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(PR2MASK)

IF("${CMAKE_SYSTEM}" MATCHES "Linux")
   set (CMAKE_CXX_STANDARD 17)
   SET(CMAKE_CXX_FLAGS "-std=c++17 -I/usr/local/include -I/usr/include")
ELSEIF("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "arm64")
   message(STATUS CMAKE_HOST_SYSTEM_PROCESSOR = ${CMAKE_HOST_SYSTEM_PROCESSOR})
   set (CMAKE_CXX_STANDARD 17)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX13.3.sdk -framework CoreFoundation")
ELSE()
   SET(CMAKE_CXX_FLAGS "-std=c++17 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -I/usr/local/include")
ENDIF()

IF(UNIX AND NOT APPLE)
  set(LINUX TRUE)
  message(STATUS ">>> Linux")
ELSE()
  message(STATUS ">>> MacOS")
ENDIF()


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

# Freetype
find_package(Freetype REQUIRED)

# target_link_libraries(your_target ${Boost_LIBRARIES})
add_executable(pr2mask pr2mask.cxx report.h report.cxx itkPolylineMask2DScanlineImageFilter.hxx)
target_include_directories(pr2mask PUBLIC ${FREETYPE_INCLUDE_DIRS})
target_link_libraries(pr2mask ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb ${FREETYPE_LIBRARIES})
#target_compile_features(pr2mask PRIVATE cxx_std_17)

add_executable(mask2rtstruct mask2rtstruct.cxx)
target_link_libraries(mask2rtstruct ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb)
#target_compile_features(mask2rtstruct PRIVATE cxx_std_17)

# imageAndMask2Report.cxx
add_executable(imageAndMask2Report imageAndMask2Report.cxx report.h report.cxx)
target_include_directories(imageAndMask2Report PUBLIC ${FREETYPE_INCLUDE_DIRS})
target_link_libraries(imageAndMask2Report ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb ${FREETYPE_LIBRARIES})
#target_compile_features(imageAndMask2Report PRIVATE cxx_std_17)

