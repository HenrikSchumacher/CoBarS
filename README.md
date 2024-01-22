# CoBarS - Conformal Barycenter Sampling

by Jason Cantarella and Henrik Schumacher

A header-only C++ library for Monte-Carlo sampling of cylic polygons with prescribed edge lengths.

# Installation

Please clone with

    git clone --recurse-submodules git@github.com:HenrikSchumacher/CoBarS.git

to load also all submodules. If you forgot to do that, you can also run the following afterwards:

    git submodule update --init --recursive
    

Pull changes from the remote repositories of any submodule by executing

    git submodule update --remote --recursive
    
# Usage

Just include the header CoBarS.hpp via

    #import "CoBarS.hpp"    

See also the examples programs in the directories Example_RandomClosedPolygon and Example_Sample_Binned for usage examples.

See also [CoBarSLink](https://github.com/HenrikSchumacher/CoBarSLink) for a more user-friendly _Mathematica_ package.
