cmake_minimum_required(VERSION 3.10)
project(chemical_visualizer)

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
if(MSVC)
    if(!EMSCRIPTEN)
	    set(_GLFW_WIN32, true)
    endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
endif(MSVC)

set(OpenGL_GL_PREFERENCE GLVND)
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

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default debug and optimize flags.
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(SOURCES_SHARED 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/assets.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/controls.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/entities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/loader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/shader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
) 


if(EMSCRIPTEN)
    add_definitions(-DEMSCRIPTEN)

    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --bind --memory-init-file 0 -s TOTAL_MEMORY=64MB")
    # Embed files at build-dir/data
    # @todo Need to figure out some way to load from file server system
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --embed-file \"${CMAKE_BINARY_DIR}/data@data\"")
    
    set(SOURCES_WASM_APP 
        ${SOURCES_SHARED}
    )    

    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/include/
    )

    add_executable(wasm_app ${SOURCES_WASM_APP})

    target_include_directories(wasm_app PUBLIC
	    ${CMAKE_CURRENT_SOURCE_DIR}/external/glm-0.9.7.1/
	    ${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include/
    )

    # Use webgl2 and glfw3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_WEBGL2=1 -s USE_GLFW=3 -s ALLOW_MEMORY_GROWTH=1")
    
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s BINARYEN_ASYNC_COMPILATION=1 -s WASM=1 -s BINARYEN=1 -s \"BINARYEN_METHOD='native-wasm'\" -s DISABLE_DEPRECATED_FIND_EVENT_TARGET_BEHAVIOR=1 -o app_wasm2.html")
    # Enable js output (no html).
    #set(CMAKE_EXECUTABLE_SUFFIX "_wasm.js")
    set(CMAKE_EXECUTABLE_SUFFIX "_emsdk.html")
    add_definitions(-DUSE_WASM)

    # Copy the content of the sourcedir/web_application directory
	add_custom_command(TARGET wasm_app POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/web_application ${CMAKE_BINARY_DIR}/
        COMMENT "Copying files from ${CMAKE_SOURCE_DIR}/web_application/ to ${CMAKE_BINARY_DIR}/")

else() # !EMSCRIPTEN

    add_subdirectory (external)

    include_directories(
	    ${CMAKE_CURRENT_SOURCE_DIR}/external/glfw-3.3.7/include/
	    ${CMAKE_CURRENT_SOURCE_DIR}/external/glm-0.9.7.1/
	    ${CMAKE_CURRENT_SOURCE_DIR}/external/glew-2.1.0/include/
	    ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/include/
	    ${CMAKE_CURRENT_SOURCE_DIR}/include/
    )

    set(ALL_LIBS
	    ${OPENGL_LIBRARY}
	    glfw
	    GLEW_210
    )

    add_definitions(
	    -DTW_STATIC
	    -DTW_NO_LIB_PRAGMA
	    -DTW_NO_DIRECT3D
	    -DGLEW_STATIC
	    -D_CRT_SECURE_NO_WARNINGS
    )

    set(SOURCES_GLFW_APP 
        ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/imgui_draw.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/imgui_impl_glfw.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/imgui_impl_opengl3.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/imgui_widgets.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/imgui_tables.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/imgui.cpp

        ${SOURCES_SHARED}
    )

    add_executable(glfw_app
        ${SOURCES_GLFW_APP}
     )

    target_link_libraries(glfw_app
	    ${ALL_LIBS}
    )

    # Copy assets to build-dir/data/ directory
    # add_custom_command(TARGET glfw_app POST_BUILD
    #     COMMAND ${CMAKE_COMMAND} -E copy_directory
    #     ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data
    #     COMMENT "Copying files from ${CMAKE_SOURCE_DIR}/data/ to ${CMAKE_BINARY_DIR}/data")

endif()