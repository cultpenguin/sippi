# GA-AEM forward
This folder contains a SIPP forward wrapper that allow calling the GA-AEM EM forward modeling code using SIPPI.

## Install ga-aem


    git clone git@github.com:GeoscienceAustralia/ga-aem.git
    cd ga-aem
    git co develop # Optional

    # On linux you will need to compile GA-AEM before use
    sudo apt-get install build-essential
    sudo apt-get install libfftw3-dev
    sudo apt-get install libopenmpi-dev

    cd ga-aem/makefiles
    export cxx=g++
    export mpicxx=mpiCC
    export cxxflags='-std=c++11 -O3 -Wall'
    export exedir='../bin/'


### Compile issues 
Please checkout https://github.com/GeoscienceAustralia/ga-aem for details on installation.


