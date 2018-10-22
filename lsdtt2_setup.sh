#!/bin/bash

# This is a startup script for LSDTopoTools
# It clones the LSDTopoTools2 repository into your home directory
# it then builds the code from there.
# Author: SMM
# Date: 15/10/2018

# Set up the inital directory
BASE_DIR=$HOME
SRC_DIR="$HOME/LSDTopoTools/LSDTopoTools2/src/"
LSD_DIR="$HOME/LSDTopoTools"
WRK_DIR="$HOME/LSDTopoTools/LSDTopoTools2"

echo "Source dir is: $SRC_DIR"
echo "LSD dir is: $LSD_DIR"

# Make the LSDTopoTools directory if it doesn't exist
if [ -f $LSD_DIR ]
  then
    echo "LSDTopoTools directory exists!"
  else
    echo "LSDTopoTools directory doesnt' exist. I'm making one"
    mkdir $LSD_DIR
fi


# clone or pull the repo, depending on what is in there
# check if the files have been cloned
if [ -f $SRC_DIR/LSDRaster.cpp ]
  then
    echo "The LSDTopoTools2 repository exists, updating to the latest version."
    git --work-tree=$WRK_DIR --git-dir=$WRK_DIR.git  pull origin master
  else
    echo "Cloning the LSDTopoTools2 repository"
    git clone https://github.com/LSDtopotools/LSDTopoTools2.git $WRK_DIR
fi

# Change the working directory to that of LSDTopoTools2/src
#echo "I am going to try to build LSDTopoTools2."
cd $SRC_DIR
echo "The current directory is:"
echo $PWD
echo "Calling the build script."
sh build.sh

# Now update the path
echo "Now I'll add the LSDTopoTools command line programs to your path."
export PATH=$HOME/LSDTopoTools/LSDTopoTools2/bin:$PATH
echo "Your path is now:"
echo $PATH
exec /bin/bash
