# GA-AEM forward
This folder contains a SIPPI forward wrapper that allow calling the GA-AEM EM forward modeling code (https://github.com/GeoscienceAustralia/ga-aem) using SIPPI.

## Download ga-aem from Github:


    git clone git@github.com:GeoscienceAustralia/ga-aem.git
  

 # On linux you will need to compile GA-AEM before use
    sudo apt-get install build-essential
    sudo apt-get install libfftw3-dev
    sudo apt-get install libopenmpi-dev

 # On linux you will need to compile GA-AEM before use
    cd ga-aem/makefiles
    export cxx=g++
    export mpicxx=mpiCC
    export cxxflags='-std=c++11 -O3 -Wall'
    export exedir='../bin/'
   
Update line 15 in gatdaem1d_matlab.make from
    
    libs       = -L$(FFTW_DIR) -lfftw3
to    

    libs       = -L$(FFTW_DIR) -lfftw -lfftw3

Then compile using 

    ./run_make.sh

### WIndows 
On windows please add libfftw3-3.dll to ga-aem/matlab/bin/x64

    cd ga-aem
    cp third_party/fftw3.2.2.dlls/64bit/libfftw3-3.dll matlab/bin/x64/.


### Compile issues 
Please checkout https://github.com/GeoscienceAustralia/ga-aem for details on installation.


