### GLFW ###
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory (glfw-3.3.7)

### OpenGl ###

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(OPENGL_LIBRARY
        ${OPENGL_LIBRARY}
        -lGL -lGLU -lXrandr -lXext -lX11 -lrt
        ${CMAKE_DL_LIBS}
        ${GLFW_LIBRARIES}
    )
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(OPENGL_LIBRARY
        ${OPENGL_LIBRARY}
        ${CMAKE_DL_LIBS}
        ${GLFW_LIBRARIES}
    )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    
        
include_directories(
    glfw-3.3.7/include/GLFW/
    glew-2.1.0/include/
)

### GLEW ###
set(GLEW_SOURCE
    glew-2.1.0/src/glew.c
)

set(GLEW_HEADERS
)

add_library( GLEW_210 STATIC
    ${GLEW_SOURCE}
    ${GLEW_INCLUDE}
)

target_link_libraries(GLEW_210
    ${OPENGL_LIBRARY}
    ${EXTRA_LIBS}
)


