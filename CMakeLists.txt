cmake_minimum_required(VERSION 3.10)
project(chemical_visualizer)

# "Pure" C++17 and no extensions.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default debug and optimize flags.
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

include_directories(${CMAKE_SOURCE_DIR})

if(NOT EMSCRIPTEN)
    # GLFW
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/subrepos/glfw/)
    find_package(OpenGL REQUIRED)
endif()


if(EMSCRIPTEN)
    add_definitions(-DEMSCRIPTEN)

    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --bind --memory-init-file 0 -s TOTAL_MEMORY=64MB")
    # Embed files at build-dir/data
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --embed-file \"${CMAKE_BINARY_DIR}/data@data\"")
    
    set(SOURCES_WASM_APP 
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/assets.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/controls.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/entities.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/graphics.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/loader.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/shader.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/utilities.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/main.cpp

    )    
    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/include/
    )

    add_executable(wasm_app ${SOURCES_WASM_APP})

    target_include_directories(wasm_app PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/subrepos/glm/
        ${CMAKE_CURRENT_SOURCE_DIR}/glad/include
    )

    # Use webgl2 and glfw3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_WEBGL2=1 -s USE_GLFW=3 -O1 -s ALLOW_MEMORY_GROWTH=1")
    
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

    # Copy assets to build-dir/data/ directory
    add_custom_command(TARGET wasm_app PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data
        COMMENT "Copying files from ${CMAKE_SOURCE_DIR}/data/ to ${CMAKE_BINARY_DIR}/data")

else()

    set(SOURCES_WASM_APP 
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/src/main.cpp
    )    
    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/glfw_app/include/
    )
    add_executable(glfw_app
        ${SOURCES_WASM_APP}
        glad/src/glad.c
    )
    target_include_directories(glfw_app PUBLIC
        ${OPENGL_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/glad/include
        ${CMAKE_CURRENT_SOURCE_DIR}/subrepos/glfw/include
        ${CMAKE_CURRENT_SOURCE_DIR}/subrepos/glm/
        ${CMAKE_CURRENT_SOURCE_DIR}/
    )
    target_link_libraries(glfw_app
        glfw
        ${GALAGA_MACHINE_LIBRARIES}
        ${OPENGL_LIBRARIES}
        ${OPENGL_gl_LIBRARY}
    )

    # Copy assets to build-dir/data/ directory
    add_custom_command(TARGET glfw_app POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data
        COMMENT "Copying files from ${CMAKE_SOURCE_DIR}/data/ to ${CMAKE_BINARY_DIR}/data")

endif()