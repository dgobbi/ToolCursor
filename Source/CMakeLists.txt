# Include directories
include_directories(${TOOLCURSOR_INCLUDE_DIRS})

# Create the main library
set(LIB_NAME vtkToolCursor)

# Sources in the current directory (library sources only!)
set(LIB_SRCS
  vtkActionCursorShapes.cxx
  vtkCursorShapes.cxx
  vtkFiducialPointsTool.cxx
  vtkFocalPlaneTool.cxx
  vtkFollowerPlane.cxx
  vtkGeometricCursorShapes.cxx
  vtkImageTool.cxx
  vtkImageToROIContourData.cxx
  vtkLassoImageTool.cxx
  vtkOpacityTool.cxx
  vtkPanCameraTool.cxx
  vtkPushPlaneTool.cxx
  vtkResliceMath.cxx
  vtkROIContourData.cxx
  vtkROIContourDataToPolyData.cxx
  vtkRotateCameraTool.cxx
  vtkSliceImageTool.cxx
  vtkSpinCameraTool.cxx
  vtkSystemCursorShapes.cxx
  vtkToolCursor.cxx
  vtkToolCursorInteractorObserver.cxx
  vtkTool.cxx
  vtkWindowLevelTool.cxx
  vtkZoomCameraTool.cxx
)

# Headers without a matching .cxx file are listed here
set(LIB_HDRS
  ${CMAKE_CURRENT_BINARY_DIR}/vtkToolCursorConfig.h
  ${CMAKE_CURRENT_BINARY_DIR}/vtkToolCursorModule.h
)

# Internal headers (not installed)
set(LIB_PRIVATE_HDRS
)

# Sources that are abstract
set(LIB_ABSTRACT
)

# Set source file properties
set_source_files_properties(${LIB_HDRS}
  PROPERTIES WRAP_EXCLUDE ON)
set_source_files_properties(${LIB_ABSTRACT}
  PROPERTIES ABSTRACT ON)
set_source_files_properties(${LIB_PRIVATE_HDRS}
  PROPERTIES SKIP_HEADER_INSTALL ON)

# Create the library
if(DEFINED VTK_MODULE_ENABLE_VTK_ToolCursor)

  # When building as a module for VTK 8.90 or later
  set(LIB_SRC_HDRS)
  foreach(_src ${LIB_SRCS})
    get_filename_component(_base ${_src} NAME_WE)
    list(APPEND LIB_SRC_HDRS ${_base}.h)
  endforeach()
  vtk_module_add_module(VTK::ToolCursor
    SOURCES ${LIB_SRCS}
    HEADERS ${LIB_SRC_HDRS} ${LIB_HDRS}
    PRIVATE_HEADERS ${LIB_PRIVATE_HDRS})

else()

  # The remainder of the file is for old, non-modular builds
  set(VTK_LIBS
    vtkFiltersHybrid
    vtkImagingCore
    vtkImagingGeneral
    vtkImagingStatistics
    vtkRenderingCore
    vtkRenderingVolume)

  set(TOOLCURSOR_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBRARIES})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vtkToolCursorModule.h.in
    "${CMAKE_CURRENT_BINARY_DIR}/vtkToolCursorModule.h" @ONLY)

  include("${TOOLCURSOR_CMAKE_DIR}/ToolCursorCommonBlock.cmake")

endif()
