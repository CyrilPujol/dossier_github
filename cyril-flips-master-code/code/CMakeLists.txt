cmake_minimum_required(VERSION 3.1...3.15)
project(demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

include_directories(BEFORE "../../include")

# .h and .cpp files
file(GLOB RESOURCE_FILES "headers/*.h" "sources/*.cpp")

# create_single_source_cgal_program("sources/demo.cpp")
add_executable(demo ${RESOURCE_FILES} ${UIS})

if(CGAL_Qt5_FOUND)
  message( STATUS "NB: CGAL_Qt5_fund !")
else()
  message( STATUS "NB: CGAL_Qt5 not fund !" )
endif()



# add headers folder(s) to the include directories
  target_include_directories(demo PUBLIC headers)

  target_link_libraries(demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)
  if(TARGET CGAL::CGAL_Core)
    target_link_libraries(demo PRIVATE CGAL::CGAL_Core)
  endif()
  if(LEDA_FOUND)
    target_link_libraries(demo PRIVATE ${LEDA_LIBRARIES})
  endif()
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS demo)
