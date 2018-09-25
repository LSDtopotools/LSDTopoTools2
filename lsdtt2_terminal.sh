# This script spawns a new shell that has the correct path to the lsdtopotools2 driver functions
# Written by Simon M Mudd
# 25-Sept-2018

# This sorts out the directory into which the code will be built
bindir=/bin
fullbindir="$PWD$bindir"
echo "The full directory for binary files is:"
echo $fullbindir

# This is the name of one of the files I will check
lsdttfile=/lsdtt-basic-metrics
fullfile="$fullbindir$lsdttfile"
echo "The name of the driver file is"
echo $fullfile

# I am checking to see if you have the binary files
if [ ! -f $fullfile ]; then
    echo "I did not find the binary files! I need to compile."
    echo "Go into the /src directory and run build.sh"
else
    echo "The binary files exist. I am updating the path and spawning a new terminal."
    echo "You can now use the lsdtt tools from this terminal."
    echo "You will need to use this script every time you want to start an LSDTopoTools session."
    export PATH=$fullbindir:$PATH   
    exec /bin/bash
fi