#!/bin/sh

# First run make to create version number form sippi.xml into version.xml
make validate


# Figure out the version number
file="version.xml" #the file where you keep your string name
if [ ! -f $file ]; then
    echo "$0: $file does not exist!"

    if [ $# -eq 0 ]
      then
        echo "$0: Set version number of first argument or run 'make' first."
      return
    else
      VERSION=$1;
    fi
else
    VERSION=$(cat "$file")

#  echo $VERSION
fi

FILENAME=SIPPI_${VERSION};
echo "$0: Creating release $FILENAME"


rm -fr SIPPI
git clone --depth 1 https://github.com/cultpenguin/sippi.git SIPPI
cd SIPPI/toolboxes
git clone --depth 1 https://github.com/cultpenguin/mgstat.git mGstat
git clone --depth 1 https://github.com/ergosimulation/mpslib.git mpslib
cd mpslib
make all
make cleano
cd ..
cd ../..

make html
mv htmldoc SIPPI/doc/.
make dblatex_ubuntu
mv sippi.pdf SIPPI/doc/.


zip -r ${FILENAME}.zip SIPPI 
tar cfz ${FILENAME}.tar.gz SIPPI
