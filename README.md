# Chemical Visualizer
View PDB (Protein Data Bank) and MOL files. This app is also built to webassembly on the [emscripten branch](https://github.com/patricklbell/chemical_visualizer/tree/emscripten).
For a brief explanation of how it works check out the acommpanying [blog post](https://patricklbell.xyz/posts/visualising-proteins).

## Haemoglobin ([PDB](https://www.rcsb.org/structure/4n7n))
![Haemoglobin Ribbon Diagram Chain Coloring](https://github.com/patricklbell/chemical_visualizer/blob/main/data/screenshots/haemoglobin_chains.png?raw=true)


## Immunoglobulin Antibody (IgG) ([PDB](https://www.rcsb.org/structure/1igt))
![Immunoglobulin Atom Diagram](https://github.com/patricklbell/chemical_visualizer/blob/main/data/screenshots/igg_atoms.png?raw=true)

## Caffeine ([MOL](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:27732))
![Caffeine Molecule Diagram](https://github.com/patricklbell/chemical_visualizer/blob/main/data/screenshots/caffeine.png?raw=true)

## Building
### Windows
#### Requirements
You need CMake and the Visual Studio build tool-chain for C++. Additionally, OpenGL
must be above version 3.3 (Should already work, if not, update your graphics 
drivers). 
#### Building
Download the source and unzip or 
```
git clone https://github.com/patricklbell/chemical_visualizer.git
```
Navigate to the root of source, create a build directory and build with CMake:
```
cd chemical_visualizer
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ./..
cmake --build . --config Release
```
CMake should copy the data folder to build/Release directory, if it doesn't copy
this yourself.

### Linux (X11)
#### Requirements
To compile GLFW you need the X11 development packages installed, on Debian and 
derivates like Ubuntu and Linux Mint the xorg-dev meta-package pulls in the 
development packages for all of X11. For more information see 
https://www.glfw.org/docs/3.3/compile.html. You will need CMake, and OpenGL 
must be above version 3.3 (Should already work, if not, update your graphics 
drivers).
#### Building
Download the source and unzip or 
```
git clone https://github.com/patricklbell/chemical_visualizer.git
```
Navigate to the root of source, create a build directory and build with CMake:
```
cd chemical_visualizer
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ./..
cmake --build .
```
CMake should copy the data folder to your build directory, if it doesn't copy 
this yourself.
```
cp -R ../data .
```

## Todo:
- MOL files which don't specify each atom's 3D coordinates but do specify bond 
angles are not displayed correctly
- Display wireframe and space-filling visualizations of PDB files
- Create different color modes for PDB files
