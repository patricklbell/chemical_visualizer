# Chemical Visualizer
View PDB (Protein Data Bank) and MOL files. This app is built with CMake both 
natively and to the web. The web build has an accompanying [demo](https://patricklbell.github.io/chemical_visualizer/). 
For a brief explanation of how the protein diagrams are generated check out [this blog post](https://patricklbell.xyz/posts/visualising-proteins).

## Haemoglobin ([PDB](https://www.rcsb.org/structure/4n7n))
![Haemoglobin Ribbon Diagram Chain Coloring](https://github.com/patricklbell/chemical_visualizer/blob/emscripten/docs/screenshots/haemoglobin_chains.png?raw=true)

## Immunoglobulin Antibody (IgG) ([PDB](https://www.rcsb.org/structure/1igt))
![Immunoglobulin Atom Diagram](https://github.com/patricklbell/chemical_visualizer/blob/emscripten/docs/screenshots/igg_atoms.png?raw=true)

## Caffeine ([MOL](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:27732))
![Caffeine Molecule Diagram](https://github.com/patricklbell/chemical_visualizer/blob/emscripten/docs/screenshots/caffeine.png?raw=true)

# Building
Download the source and unzip or 
```bash
    git clone https://github.com/patricklbell/chemical_visualizer.git
```

## Web
The build system uses CMake and Emscripten to create a WebAssembly module, 
it is easiest to use inside the accompanying Docker container. For native 
building see instructions below.

### With VSCode using Docker

The provided [Dockerfile](.devcontainer/Dockerfile) contains all you need to 
build and run this project (e.g. Emscripten SDK version 2.0+, gcc, CMake, Ninja).


Open the folder with VSCode using "Remote-Container: Open folder in container", 
Enable the CMake Tools extension for VSCode.

### Building the WebAssembly based application with Emscripten SDK

Create a running container from the [Dockerfile](.devcontainer/Dockerfile), mount 
the root of the repository and navigate into the workspace, and run:

  ```bash
    mkdir build && cd build && ln -s ../data data
    emcmake cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
    cmake --build .
  ```

You can replace Ninja with any project file generator you like, e.g. make/VisualStudio/XCode/Eclipse.

## Natively
### Windows
#### Requirements
You need CMake and the Visual Studio build tool-chain for C++. Additionally, OpenGL
must be above version 3.3 (Should already work, if not, update your graphics 
drivers). 
#### Build
Open the root of source with Visual Studio as a folder, right click CMakeList.txt and select 'Set as Startup Item.' 
Set the CWD to the root of the project, or copy /data to the build folder. Build

### Linux (X11)
#### Requirements
To compile GLFW you need the X11 development packages installed, on Debian and 
derivates like Ubuntu and Linux Mint the xorg-dev meta-package pulls in the 
development packages for all of X11. For more information see 
https://www.glfw.org/docs/3.3/compile.html. You will need CMake, and OpenGL 
must be above version 3.3 (Should already work, if not, update your graphics 
drivers).
#### Build
Navigate to the root of source, create a build directory, link data, and build with CMake:
```bash
    cd chemical_visualizer
    mkdir build && cd build && ln -s ../data data
    cmake -DCMAKE_BUILD_TYPE=Release ..
    cmake --build .
```