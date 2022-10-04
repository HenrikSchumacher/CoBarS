# CyclicSampler

by Jason Cantarella and Henrik Schumacher

A program for Monte-Carlo sampling of cylic polygons with prescribed edge lengths.

# Installation

After cloning make sure to run

    git submodule init update --recursive
    

Pull changes from the remote repositories of any submodule by executing

    git submodule update --remote --recursive
    
    
If you accidentally modified one of the submodules you can run

    git submodule foreach --recursive git reset --hard
    
to repair this.
