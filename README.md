# Chemical Visualizer
View PDB (Protein Data Bank) and mol files. This is the emscripten branch, which builts with cmake and emcc to create a static website.
See the result at https://patricklbell.github.io/chemical_visualizer/.

## Building
The build system is modified from https://github.com/lukka/CppOpenGLWebAssemblyCMake and requires Docker. Note that while the build system supports building natively, this is untested (won't work) since main already builds natively.

### With VSCode using Docker

The provided [Dockerfile](.devcontainer/Dockerfile) contains all you need to build and run this project (e.g. Emscripten SDK version 2.0+, gcc, CMake, Ninja).

Open the folder with VSCode using "Remote-Container: Open folder in container", 

Enable the CMake Tools extension for VSCode.

### Building the WebAssembly based application with Emscripten SDK

Create a running container from the [Dockerfile](.devcontainer/Dockerfile), mount the root of the repository onto /workspace/chemical_visualizer_emscripten/,
and run:

  ```bash
    cd /workspace/chemical_visualizer_emscripten/
    mkdir build && cd build
    emcmake cmake -GNinja ..
    cmake --build .
  ```

### Building the native application for Linux/macOS/Windows

  ```bash
    mkdir build && cd build
    cmake -GNinja ..
    cmake --build .
  ```

You can replace Ninja with any project file generator you like, e.g. make/VisualStudio/XCode/Eclipse.