# CycleSampler

by Jason Cantarella and Henrik Schumacher

A program for Monte-Carlo sampling of cylic polygons with prescribed edge lengths.

# Installation

Please clone with

    git clone --recurse-submodules git@github.com:HenrikSchumacher/CycleSampler.git

to load also all submodules. Id you forgot to do that, you can also run the following afterwards:

    git submodule update --init --recursive
    

Pull changes from the remote repositories of any submodule by executing

    git submodule update --remote --recursive
    
    
# Trouble shooting

If you accidentally modified one of the submodules you can run

    git submodule foreach --recursive git reset --hard
    
to repair this.
