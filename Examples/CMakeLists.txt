project(TOOLCURSORExamples)

# Building outside of TOOLCURSOR build tree
if(NOT TOOLCURSOR_BINARY_DIR)
  message(STATUS "Buiding examples outside of TOOLCURSOR build tree")

  IF(NOT TOOLCURSOR_DIR)
    find_package(TOOLCURSOR REQUIRED)
  else(NOT TOOLCURSOR_DIR)
    include(${TOOLCURSOR_DIR}/TOOLCURSORConfig.cmake)
  endif()

  include(${TOOLCURSOR_USE_FILE})

  # Test output
  set(TOOLCURSOR_EXAMPLES_TEMP_DIR ${TOOLCURSOR_EXAMPLES_BINARY_DIR}/Testing/Temporary)

  # Output paths
  set(EXECUTABLE_OUTPUT_PATH ${TOOLCURSOR_EXAMPLES_BINARY_DIR}/bin/ CACHE INTERNAL "Single output directory for building all executables.")
  set(TOOLCURSOR_EXAMPLES_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
  set(CXX_TEST_PATH ${TOOLCURSOR_EXAMPLES_EXECUTABLE_DIR})

# Testing inside the TOOLCURSOR build tree
else(NOT TOOLCURSOR_BINARY_DIR)
  message(STATUS "Buiding examples inside the ToolCursor build tree")

  include_directories(${TOOLCURSOR_INCLUDE_DIRS})

  # Test output
  set(TOOLCURSOR_EXAMPLES_TEMP_DIR ${TOOLCURSOR_BINARY_DIR}/Testing/Temporary)
endif()

# The library dependencies
set(VTK_TARGET_PREFIX "VTK::")
if(VTK_VERSION VERSION_LESS 8.90)
  set(VTK_TARGET_PREFIX "vtk")
endif()
set(VTK_LIBS)

# required libraries
foreach(_lib RenderingCore RenderingImage InteractionStyle
  RenderingFreeType IOImage IOMINC CommonSystem CommonCore
  ImagingStencil ImagingMath FiltersSources FiltersGeneral)

  list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
endforeach()

# factory libraries (only needed if present)
foreach(_lib IOMPIImage
  RenderingOpenGL RenderingFreeTypeOpenGL RenderingVolumeOpenGL
  RenderingOpenGL2 RenderingFreeTypeOpenGL2 RenderingVolumeOpenGL2)

  list(FIND VTK_LIBRARIES ${VTK_TARGET_PREFIX}${_lib} _index)
  if(_index GREATER -1)
    list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
  endif()
endforeach()

# utility libraries
if(VTK_VERSION VERSION_LESS 8.90)
  list(APPEND VTK_LIBS vtksys)
else()
  list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}vtksys)
endif()

add_executable(LassoImageTool LassoImageTool.cxx)
target_link_libraries(LassoImageTool
  ${VTK_TARGET_PREFIX}ToolCursor ${VTK_LIBS})

if(NOT VTK_VERSION VERSION_LESS 8.90)
  vtk_module_autoinit(
    TARGETS LassoImageTool
    MODULES VTK::RenderingOpenGL2 VTK::RenderingUI VTK::InteractionStyle)
endif()
