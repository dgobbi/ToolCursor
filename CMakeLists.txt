cmake_minimum_required(VERSION 3.12)

project(TOOLCURSOR)

# Allow custom cmake overrides
set(TOOLCURSOR_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
set(CMAKE_MODULE_PATH "${TOOLCURSOR_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

# Project version
set(TOOLCURSOR_MAJOR_VERSION 1)
set(TOOLCURSOR_MINOR_VERSION 0)
set(TOOLCURSOR_PATCH_VERSION 0)
set(TOOLCURSOR_SHORT_VERSION "${TOOLCURSOR_MAJOR_VERSION}.${TOOLCURSOR_MINOR_VERSION}")
set(TOOLCURSOR_VERSION "${TOOLCURSOR_SHORT_VERSION}.${TOOLCURSOR_PATCH_VERSION}")

# Include directories
set(TOOLCURSOR_INCLUDE_DIRS)
if(VTK_VERSION VERSION_LESS 8.90)
  list(APPEND TOOLCURSOR_INCLUDE_DIRS
    "${TOOLCURSOR_SOURCE_DIR}/Source"
    "${TOOLCURSOR_BINARY_DIR}/Source")
endif()

# Find VTK and get VTK configuration information
find_package(VTK REQUIRED)

# Define VTK_VERSION if not already defined
if(NOT DEFINED VTK_VERSION)
  set(VTK_VERSION ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})
endif()
if(VTK_VERSION VERSION_LESS 8.90)
  include(${VTK_USE_FILE})
endif()
if(VTK_VERSION VERSION_LESS 7.1)
  message(FATAL_ERROR "VTK 7.1 or later is required.")
endif()

# Standard option
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

# Use ctest for dashboard submissions
option(BUILD_TESTING "Build the tests" OFF)
if(BUILD_TESTING)
  include(CTest)
endif()

# Which parts of the package to build
option(BUILD_EXAMPLES "Build the examples" OFF)

# Wrap Python option
if(VTK_WRAP_PYTHON)
  option(BUILD_PYTHON_WRAPPERS "Build python wrappers" ON)
else()
  unset(BUILD_PYTHON_WRAPPERS CACHE)
endif()

# Configuration header
set(TOOLCURSOR_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBRARIES})
CONFIGURE_FILE(CMake/vtkToolCursorConfig.h.in
  "${TOOLCURSOR_BINARY_DIR}/Source/vtkToolCursorConfig.h" @ONLY)

# Defines CMAKE_INSTALL_LIBDIR etcetera
include(GNUInstallDirs)

# Set up our directory structure for output libraries and binaries
# (Note: these are the build locations, not the install locations)
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${TOOLCURSOR_BINARY_DIR}/bin")
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  if(UNIX)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${TOOLCURSOR_BINARY_DIR}/lib")
  else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${TOOLCURSOR_BINARY_DIR}/bin")
  endif()
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${TOOLCURSOR_BINARY_DIR}/lib")
endif()

# Building for VTK 9.x or later
if(NOT VTK_VERSION VERSION_LESS 8.90)
  # Clear VTK_LIBS, linking will be handled by modular build
  set(VTK_LIBS)

  option(TOOLCURSOR_VERSIONED_INSTALL "Install with versioned names." ON)
  mark_as_advanced(TOOLCURSOR_VERSIONED_INSTALL)
  set(vtk_version_suffix "")
  set(toolcursor_library_suffix "")
  if(VTK_VERSIONED_INSTALL OR TOOLCURSOR_VERSIONED_INSTALL)
    set(vtk_version_suffix "-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
    set(toolcursor_library_suffix "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  endif()
  if(DEFINED TOOLCURSOR_CUSTOM_LIBRARY_SUFFIX)
    set(toolcursor_library_suffix "${TOOLCURSOR_CUSTOM_LIBRARY_SUFFIX}")
  endif()

  # Optional arguments for vtk_module_build()
  set(module_build_extra_args)
  if(NOT VTK_VERSION VERSION_LESS 9.4.0)
    list(APPEND module_build_extra_args
      USE_FILE_SETS ON)
  endif()
  if(NOT "${toolcursor_library_suffix}" STREQUAL "")
    list(APPEND module_build_extra_args
      LIBRARY_NAME_SUFFIX "${toolcursor_library_suffix}")
  endif()

  # Where to put config files
  set(cmake_destination
    "${CMAKE_INSTALL_LIBDIR}/cmake/toolcursor-${TOOLCURSOR_SHORT_VERSION}")

  set(toolcursor_module_files
    "${CMAKE_CURRENT_SOURCE_DIR}/Source/vtk.module")

  vtk_module_scan(
    MODULE_FILES        ${toolcursor_module_files}
    REQUEST_MODULES     "VTK::ToolCursor"
    PROVIDES_MODULES    toolcursor_modules
    ENABLE_TESTS        "${BUILD_TESTING}")

  vtk_module_build(
    MODULES             ${toolcursor_modules}
    INSTALL_EXPORT      TOOLCURSOR
    ARCHIVE_DESTINATION "lib"
    HEADERS_DESTINATION "include/vtk${vtk_version_suffix}"
    CMAKE_DESTINATION   "${cmake_destination}"
    LICENSE_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/toolcursor-${TOOLCURSOR_SHORT_VERSION}"
    HIERARCHY_DESTINATION "${CMAKE_INSTALL_LIBDIR}/vtk${vtk_version_suffix}/hierarchy/ToolCursor"
    VERSION             "${TOOLCURSOR_VERSION}"
    SOVERSION           "1"
    ${module_build_extra_args})

  # Disable deprecation warnings for MSVC
  if(MSVC)
    vtk_module_definitions(${toolcursor_modules}
      PRIVATE -D_CRT_SECURE_NO_DEPRECATE)
  endif()

  # Set targets to be exported by toolcursor-config.cmake
  set(TOOLCURSOR_LINK_TARGETS ${toolcursor_modules})

  if(BUILD_PYTHON_WRAPPERS)
    find_package(Python${VTK_PYTHON_VERSION} QUIET COMPONENTS Development)

    vtk_module_wrap_python(
      MODULES         ${toolcursor_modules}
      TARGET          VTK::toolcursormodules
      INSTALL_EXPORT  TOOLCURSORPython
      PYTHON_PACKAGE  "toolcursor"
      CMAKE_DESTINATION   "${cmake_destination}"
      LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      SOABI           "${Python${VTK_PYTHON_VERSION}_SOABI}"
      BUILD_STATIC    OFF)

    vtk_module_python_default_destination(python_destination)

    export(
      EXPORT    TOOLCURSORPython
      NAMESPACE VTK::
      FILE      "${cmake_destination}/TOOLCURSORPython-targets.cmake")
    install(
      EXPORT      TOOLCURSORPython
      NAMESPACE   VTK::
      FILE        TOOLCURSORPython-targets.cmake
      DESTINATION "${cmake_destination}")
  endif()

  # Generate config files
  include(CMakePackageConfigHelpers)
  configure_package_config_file(CMake/toolcursor-config.cmake.in
    "${cmake_destination}/toolcursor-config.cmake"
    INSTALL_DESTINATION "${cmake_destination}")
  write_basic_package_version_file("${cmake_destination}/toolcursor-config-version.cmake"
    VERSION "${TOOLCURSOR_MAJOR_VERSION}.${TOOLCURSOR_MINOR_VERSION}.${TOOLCURSOR_PATCH_VERSION}"
    COMPATIBILITY AnyNewerVersion)

  # In root of build directory for backwards compatibility
  file(GENERATE
    OUTPUT  "${CMAKE_BINARY_DIR}/toolcursor-config.cmake"
    CONTENT "include(\"${CMAKE_BINARY_DIR}/${cmake_destination}/toolcursor-config.cmake\")\n")
  configure_file(
    "${CMAKE_BINARY_DIR}/${cmake_destination}/toolcursor-config-version.cmake"
    "${CMAKE_BINARY_DIR}/toolcursor-config-version.cmake"
    COPYONLY)

  install(
    FILES       "${CMAKE_BINARY_DIR}/${cmake_destination}/toolcursor-config.cmake"
                "${CMAKE_BINARY_DIR}/${cmake_destination}/toolcursor-config-version.cmake"
    DESTINATION "${cmake_destination}"
    COMPONENT   "development")

endif() #(NOT VTK_VERSION VERSION_LESS 8.90)

# Set up our target directory structure for "make install"
set(TOOLCURSOR_BIN_DIR "bin") # for executables and ".dll" libraries
if(NOT TOOLCURSOR_LIB_DIR)
  set(TOOLCURSOR_LIB_DIR "lib") # for unix/linux/osx shared libraries
endif()
if(NOT TOOLCURSOR_ARC_DIR)
  set(TOOLCURSOR_ARC_DIR "${TOOLCURSOR_LIB_DIR}") # for static or ".lib" libraries
endif()
if(NOT TOOLCURSOR_INC_DIR)
  set(TOOLCURSOR_INC_DIR "include") # for header files
endif()
if(NOT TOOLCURSOR_DATA_DIR)
  set(TOOLCURSOR_DATA_DIR "share/toolcursor-${TOOLCURSOR_SHORT_VERSION}")
endif()
if(NOT TOOLCURSOR_PGK_DIR)
  set(TOOLCURSOR_PGK_DIR "${TOOLCURSOR_LIB_DIR}/toolcursor-${TOOLCURSOR_SHORT_VERSION}/cmake")
endif()

# Set the full paths to the install tree
set(TOOLCURSOR_RUNTIME_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${TOOLCURSOR_BIN_DIR})
set(TOOLCURSOR_LIBRARY_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${TOOLCURSOR_LIB_DIR})
set(TOOLCURSOR_ARCHIVE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${TOOLCURSOR_ARC_DIR})
set(TOOLCURSOR_INCLUDE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${TOOLCURSOR_INC_DIR})
set(TOOLCURSOR_DATA_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${TOOLCURSOR_DATA_DIR})
set(TOOLCURSOR_CMAKE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${TOOLCURSOR_PGK_DIR})

# Libraries
set(TOOLCURSOR_LIBRARIES vtkToolCursor)

if(VTK_VERSION VERSION_LESS 8.90)
  # External wrapping previous to VTK 9

  set(VTK_WRAP_INCLUDE_DIRS
    ${TOOLCURSOR_INCLUDE_DIRS}
    ${VTK_WRAP_INCLUDE_DIRS})

  # For the python wrappers
  if(BUILD_PYTHON_WRAPPERS)
    # TODO: better python package destination for VTK 7 and 8
    set(python_destination "lib")
    if(TARGET vtkCommonCorePythonD OR TARGET vtkCommonKitPythonD)
      set(TOOLCURSOR_PYTHON_LIBRARIES vtkToolCursorPythonD)
    endif()
  endif()
endif()

if(BUILD_PYTHON_WRAPPERS)
  string(CONCAT _content
    "from .vtkToolCursor import *\n\n"
    "__all__ = ['vtkToolCursor']\n"
    "__version__ = \"${TOOLCURSOR_VERSION}\"\n")
  if(NOT VTK_VERSION VERSION_LESS 8.90)
    string(PREPEND _content "import vtkmodules\n")
  endif()
  file(GENERATE
    OUTPUT  "${CMAKE_BINARY_DIR}/${python_destination}/toolcursor/__init__.py"
    CONTENT "${_content}")
  install(
    FILES       "${CMAKE_BINARY_DIR}/${python_destination}/toolcursor/__init__.py"
    DESTINATION "${python_destination}/toolcursor/")
endif()

# Disable C and C++ deprecation warnings for MSVC 2005 and later
if(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400 OR MSVC10)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE
                  -D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_SCL_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
endif()

# Add subdirs if not handled by modular build
if(VTK_VERSION VERSION_LESS 8.90)
  add_subdirectory(Source)
endif()

# Subdirectories
if(BUILD_EXAMPLES)
  add_subdirectory(Examples)
endif()

#if(BUILD_TESTING)
#  enable_testing()
#  add_subdirectory(Testing)
#endif()

# Install files
install(FILES Copyright.txt
  DESTINATION ${TOOLCURSOR_DATA_INSTALL_DEST})

if(VTK_VERSION VERSION_LESS 8.90)
  #---------------------------------------------------
  # Export information for other projects

  # Add all targets to the build-tree export set
  export(TARGETS ${TOOLCURSOR_LIBRARIES}
    ${TOOLCURSOR_PYTHON_LIBRARIES}
    FILE "${TOOLCURSOR_BINARY_DIR}/ToolCursorTargets.cmake")

  # Export the package for use from the build-tree
  # (this registers the build-tree with a global CMake-registry)
  export(PACKAGE TOOLCURSOR)

  # Create the ToolCursorConfig.cmake and ToolCursorConfigVersion.cmake files
  # ... for the build tree
  set(INCLUDE_DIRS_CONFIG ${TOOLCURSOR_INCLUDE_DIRS})
  set(LIBRARY_DIRS_CONFIG "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
  configure_file(CMake/ToolCursorConfig.cmake.in
    "${TOOLCURSOR_BINARY_DIR}/ToolCursorConfig.cmake" @ONLY)
  # ... for the install tree
  file(RELATIVE_PATH REL_INCLUDE_DIR "${TOOLCURSOR_CMAKE_INSTALL_DEST}"
    "${TOOLCURSOR_INCLUDE_INSTALL_DEST}")
  set(INCLUDE_DIRS_CONFIG "\${TOOLCURSOR_PGK_DIR}/${REL_INCLUDE_DIR}")
  file(RELATIVE_PATH REL_LIBRARY_DIR "${TOOLCURSOR_CMAKE_INSTALL_DEST}"
    "${TOOLCURSOR_LIBRARY_INSTALL_DEST}")
  set(LIBRARY_DIRS_CONFIG "\${TOOLCURSOR_PGK_DIR}/${REL_LIBRARY_DIR}")
  configure_file(CMake/ToolCursorConfig.cmake.in
    "${TOOLCURSOR_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ToolCursorConfig.cmake" @ONLY)

  # ... for both
  configure_file(CMake/ToolCursorConfigVersion.cmake.in
    "${TOOLCURSOR_BINARY_DIR}/ToolCursorConfigVersion.cmake" @ONLY)

  # Install the ToolCursorConfig.cmake and ToolCursorConfigVersion.cmake
  install(FILES
    "${TOOLCURSOR_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ToolCursorConfig.cmake"
    "${TOOLCURSOR_BINARY_DIR}/ToolCursorConfigVersion.cmake"
    DESTINATION "${TOOLCURSOR_CMAKE_INSTALL_DEST}" COMPONENT Development)

  # Install the Use file in the build and install directories
  configure_file(${TOOLCURSOR_SOURCE_DIR}/CMake/UseToolCursor.cmake.in
    ${TOOLCURSOR_BINARY_DIR}/UseToolCursor.cmake COPYONLY IMMEDIATE)
  install(FILES
    "${TOOLCURSOR_BINARY_DIR}/UseToolCursor.cmake"
    DESTINATION "${TOOLCURSOR_CMAKE_INSTALL_DEST}" COMPONENT Development)

  # Install the export set for use with the install-tree
  install(EXPORT ToolCursorTargets
    DESTINATION "${TOOLCURSOR_CMAKE_INSTALL_DEST}" COMPONENT Development)
endif()
