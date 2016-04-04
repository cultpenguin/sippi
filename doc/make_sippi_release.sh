#!/bin/sh

if [ $# -eq 0 ]
  then
    echo "$0: No arguments supplied (must be version number)"

  return	
fi

FILENAME=SIPPI_$1.zip;

rm -fr SIPPI

git clone --depth 1 https://github.com/cultpenguin/sippi.git SIPPI
cd SIPPI/toolboxes
git clone --depth 1 https://github.com/cultpenguin/mgstat.git mGstat
git clone --depth 1 https://github.com/ergosimulation/mps.git MPS
cd ../..
zip -r $FILENAME SIPPI




