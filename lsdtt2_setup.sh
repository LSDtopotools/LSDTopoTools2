#!/bin/bash

# This is a startup script for LSDTopoTools
# It clones the LSDTopoTools2 repository into a directory of your choice
# it then builds the code from there.
# Author: SMM
# Date: 22/05/2020

echo "Hello, I am going to set up LSDTopoTools2 for you."


echo "First, would you like me to see if you already have an LSDTopoTools2 directory? "
echo "This uses a find function so it might take a while."

while true; do
  echo "Enter y or n"
  read yn
  case $yn in
    [Yy]* ) echo "Checking...";
            find / -type d -name "LSDTopoTools2";
            echo "end check.";
            echo "If a directory has been found it means you might already have LSDTopoTools2."
            echo "If you don't want to install LSDTopoTools2 twice,"
            echo "use this directory as the installation path in the options below."                         
            sleep 3
            break;;
    [Nn]* ) echo "Okay, I'll move on"; break;;
    * ) echo "Please answer yes or no.";;
  esac
done

echo "===================================================================="
echo "Okay, now, I'd like to know where you want me to set up."
echo "Your current directory is: "
echo $PWD
echo "Where yould you like me to put LSDTopoTools?"
echo "I will make an LSDTopoTools directory here if it doesn't exist"
echo "1 -- I am in an LSDTopoTools docker container. Make all the choices for me."
echo "2 -- Put it in this directory."
echo "3 -- Put it in my home directory."
echo "4 -- Put it somewhere else. I'll ask you where you want it if you make this choice."
echo "5 -- I have conda. Just install it. Note that this doesn't create any of the test data directories"
read choice


echo "You chose:"
echo $choice
ROOT_DIR="/root"

if [ $choice -eq 1 ]; then
  echo "You are in an LSDTopoTools docker container. I am going to make the choices for you."
  BASE_DIR=$HOME
  echo "The base directory is"
  echo $BASE_DIR
  if [ $BASE_DIR = $ROOT_DIR ]; then  
    BASE_DIR=""
  fi
elif [ $choice -eq 2 ]; then
  echo "You chose this directory"
  BASE_DIR=$PWD
elif [ $choice -eq 3 ]; then
  echo "I will put LSDTopoTools2 in your home directory"
  BASE_DIR=$HOME
  if [ $BASE_DIR = $ROOT_DIR ]; then  
    BASE_DIR=""
  fi
elif [ $choice -eq 4 ]; then
  echo "You chose somewhere else."
  echo "Please enter the full directory path, without trailing slash, where you want "
  read choice_dir
  echo "Checking to see if"
  echo $choice_dir
  echo "exists"
  if [ -f $choice_dir ]; then
    echo "The directory exists. Setting that as the LSDTopoTools base directory"
    BASE_DIR=$choice_dir
  else
    echo "The directory doesn't exist. I am exiting. Please create your directory of choice using mkdir"
    exit 1
  fi
elif [ $choice -eq 5 ]; then
  echo "You have conda? Excellent. I will just install the tools."
  if conda ; then
    conda install -y -c conda-forge lsdtopotools
  else
    echo "You don't have conda installed. Exiting."
    exit 1
  fi
else
  echo "Sorry, I didn't understand that. Please try again."
  exit 1
fi

echo "You have chosen to install LSDTopoTools here:"
echo $BASE_DIR

# Now set up all the subsidiary directories
SRC_DIR="$BASE_DIR/LSDTopoTools/LSDTopoTools2/src/"
LSD_DIR="$BASE_DIR/LSDTopoTools"
WRK_DIR="$BASE_DIR/LSDTopoTools/LSDTopoTools2"
DATA_DIR="$BASE_DIR/LSDTopoTools/data/ExampleTopoDatasets"

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


echo "====================================================================="
echo "Would you like me to get the Example data?"
echo "It is ~200Mb"

while true; do
  echo "Enter y or n"
  read yn
  case $yn in
    [Yy]* ) echo "Okay, checking for data...";
            if [ -f $DATA_DIR ]
              then
                echo "The Example data already exists!."
              else
                echo "Cloning the LSDTopoTools2 repository"
                git clone https://github.com/LSDtopotools/ExampleTopoDatasets.git $DATA_DIR
            fi
            break;;
    [Nn]* ) echo "Okay, I'll move on"; break;;
    * ) echo "Please answer yes or no.";;
  esac
done



# Now update the path
echo "Now I'll add the LSDTopoTools command line programs to your path."
export PATH=$BASE_DIR/LSDTopoTools/LSDTopoTools2/bin:$PATH
echo "Your path is now:"
echo $PATH
exec /bin/bash


