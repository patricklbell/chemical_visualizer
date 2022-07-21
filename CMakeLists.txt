# CMake entry point
cmake_minimum_required (VERSION 3.0)
project (chem_viz)

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
if(MSVC)
	set(_GLFW_WIN32, true)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
endif(MSVC)

find_package(OpenGL 33 REQUIRED)

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
	message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
	message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()

# Compile external dependencies 
add_subdirectory (external)

include_directories(
	external/glfw-3.3.7/include/
	external/glm-0.9.7.1/
	external/glew-2.1.0/include/
	external/imgui/include/
	external/stb/
	include/
)

set(ALL_LIBS
	${OPENGL_LIBRARY}
	glfw
	GLEW_210
)

option(USE_ASSIMP "Compile with assimp to load obj not custom mesh file" OFF)

if(USE_ASSIMP)
    include_directories(
	    external/assimp-5.1.5/include/
    )
    list(APPEND ALL_LIBS assimp)
endif()

add_definitions(
	-DTW_STATIC
	-DTW_NO_LIB_PRAGMA
	-DTW_NO_DIRECT3D
	-DGLEW_STATIC
	-D_CRT_SECURE_NO_WARNINGS
)


set(SOURCES 
    external/imgui/imgui_draw.cpp
    external/imgui/imgui_impl_glfw.cpp
    external/imgui/imgui_impl_opengl3.cpp
    external/imgui/imgui_widgets.cpp
    external/imgui/imgui_tables.cpp
    external/imgui/imgui.cpp

	src/main.cpp
	src/shader.cpp
	src/texture.cpp
	src/utilities.cpp
	src/controls.cpp
	src/editor.cpp
	src/graphics.cpp
	src/assets.cpp
	src/entities.cpp
	src/loader.cpp
)

# YCM generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE Debug)

add_executable (chem_viz ${SOURCES})

target_link_libraries(chem_viz 
	${ALL_LIBS}
)

if(USE_ASSIMP)
    set_target_properties(chem_viz PROPERTIES COMPILE_DEFINITIONS "USE_ASSIMP")
endif()

add_custom_target(copy-runtime-files ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data)

# Copy compile_commands.json
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()
