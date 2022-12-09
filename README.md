# CycleSampler

by Jason Cantarella and Henrik Schumacher

A header-only C++ library for Monte-Carlo sampling of cylic polygons with prescribed edge lengths.

# Installation

Please clone with

    git clone --recurse-submodules git@github.com:HenrikSchumacher/CycleSampler.git

to load also all submodules. If you forgot to do that, you can also run the following afterwards:

    git submodule update --init --recursive
    

Pull changes from the remote repositories of any submodule by executing

    git submodule update --remote --recursive
    
Currently, the package depends on eigen (see https://eigen.tuxfamily.org) and OpenMP. So please make sure that they are installed and found by the compiler. So far it has been tested only under macos with Apple Clang as compiler. But it should compile equally fine on other architectures and with other compilers.
    
# Usage

Just include the header CycleSampler.hpp. See also the examples programs in the directories Example_RandomClosedPolygon and Example_Sample_Binned for usage examples.
    
# Trouble shooting

If you accidentally modified one of the submodules you can run

    git submodule foreach --recursive git reset --hard
    
to repair this.
