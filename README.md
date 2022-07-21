# Chemical Visualizer
View PDB (Protein Data Bank) and mol files.

# 8DH6 ([PDB](https://www.rcsb.org/structure/8DH6))
![8DH6 PDB file](https://github.com/patricklbell/chemical_visualizer/blob/main/screenshot_pdb_8dh6.png?raw=true)

# Caffeine (Molfile)
![Caffeine Molfile](https://github.com/patricklbell/chemical_visualizer/blob/main/screenshot_mol_caffeine.png?raw=true)

# Building
## Windows
### Requirements
You need CMake and the Visual Studio build toochain for C++. Additionally OpenGL must be above version 3.3 (Should already work, if not, update your graphics drivers). 
### Building
Download the source and unzip or 
```
git clone https://github.com/patricklbell/chemical_visualizer.git
```
Navigate to the root of source, create a build directory and build with CMake:
```
mkdir build
cd build
cmake ./..
cmake --build .
```
CMake should copy the data folder to your build directory, if it doesn't copy this yourself.
```
cp -R ../data .
```
## Linux (X11)
### Requirements
To compile GLFW you need the X11 development packages installed, on Debian and derivates like Ubuntu and Linux Mint the xorg-dev meta-package pulls in the development packages for all of X11. For more information see https://www.glfw.org/docs/3.3/compile.html. You will need CMake and OpenGL must be above version 3.3 (Should already work, if not, update your graphics drivers).
### Building
Download the source and unzip or 
```
git clone https://github.com/patricklbell/chemical_visualizer.git
```
Navigate to the root of source, create a build directory and build with CMake:
```
mkdir build
cd build
cmake ./..
cmake --build .
```
CMake should copy the data folder to your build directory, if it doesn't copy this yourself.
```
cp -R ../data .
```