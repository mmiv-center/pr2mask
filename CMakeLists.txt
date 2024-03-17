cmake_minimum_required(VERSION 3.20.1)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(PR2MASK)

if (EXISTS "/Users/hauke/src/itk/InsightToolkit-5.3.0/bin")
  set(ITK_DIR /Users/hauke/src/itk/InsightToolkit-5.3.0/bin)
else()
  set(ITK_DIR /opt/InsightToolkit-5.3.0/bin)
endif()

#set(CMAKE_CXX_STANDARD 20)
#set(Qt5_DIR /usr/local/Cellar/qt@5/5.15.10/lib/cmake/Qt5)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -framework CoreFoundation")

# find the correct path for isysroot with `xcrun --show-sdk-path`
IF("${CMAKE_SYSTEM}" MATCHES "Linux")
   if (EXISTS "/opt/InsightToolkit-5.3.0/bin")
      set(ITK_DIR /opt/InsightToolkit-5.3.0/bin)
   endif()

   set (CMAKE_CXX_STANDARD 17)
   SET(CMAKE_CXX_FLAGS "-std=c++17 -I/usr/local/include -I/usr/include")
ELSEIF("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "arm64")
   # InsightToolkit-5.3.0			
   if (EXISTS "../InsightToolkit-5.3.0/bin")
      set(ITK_DIR ..//InsightToolkit-5.3.0/bin)
   endif()
   message(STATUS CMAKE_HOST_SYSTEM_PROCESSOR = ${CMAKE_HOST_SYSTEM_PROCESSOR})
   set (CMAKE_CXX_STANDARD 17)
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
set(boost_libs system filesystem timer date_time program_options)
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

set( GDCM_LIBRARIES
    gdcmCommon gdcmDICT gdcmDSED gdcmIOD gdcmMEXD gdcmMSFF gdcmcharls
    gdcmjpeg12 gdcmjpeg16 gdcmjpeg8 ITKIOGDCM
)

if (EXISTS "../dcmtk")
  set(DCMTK_DIR ../dcmtk)
  MESSAGE(STATUS ">>> Found local DCMTK")
endif()

find_package(DCMTK REQUIRED)

find_package(LIBZIP REQUIRED)
if (LIBZIP_FOUND)
   MESSAGE(STATUS ">>> Found LIBZIP")
ELSE(LIBZIP_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build app without libzip.")
ENDIF()


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
target_link_libraries(imageAndMask2Report ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb ${FREETYPE_LIBRARIES} ${GDCM_LIBRARIES} libzip::zip)
target_compile_features(imageAndMask2Report PRIVATE cxx_std_17)

# structured report from json
add_executable(json2SR json2SR.cxx)
target_include_directories(json2SR PUBLIC ${DCMTK_INCLUDE_DIRS})
target_link_directories(json2SR PUBLIC ../dcmtk/dcmsr/libsrc ../dcmtk/dcmimage/libsrc ../dcmtk/dcmimgle/libsrc ../dcmtk/dcmdata/libsrc ../dcmtk/oflog/libsrc ../dcmtk/ofstd/libsrc ../dcmtk/oficonv/libsrc ../dcmtk/lib)
target_link_libraries(json2SR dcmsr dcmimage dcmimgle dcmdata oflog ofstd oficonv xml2 z zlib ${Boost_LIBRARIES} pthread)
#target_compile_features(json2SR PRIVATE cxx_std_17)

# Interpolate labels in a mask file with gaps
add_executable(MorphologicalContourInterpolation MorphologicalContourInterpolation.cxx)
#target_include_directories(MorphologicalContourInterpolation PUBLIC ${DCMTK_INCLUDE_DIRS})
#target_link_directories(MorphologicalContourInterpolation PUBLIC ../dcmtk/dcmsr/libsrc ../dcmtk/dcmimage/libsrc ../dcmtk/dcmimgle/libsrc ../dcmtk/dcmdata/libsrc ../dcmtk/oflog/libsrc ../dcmtk/ofstd/libsrc ../dcmtk/oficonv/libsrc ../dcmtk/lib)
target_link_libraries(MorphologicalContourInterpolation ${ITK_LIBRARIES} ${Boost_LIBRARIES} pthread)
