export LDFLAGS="-Wl,-rpath,$PREFIX/lib -L$PREFIX/lib $LDFLAGS"

# This sorts out the directory into which the code will be built
dir=$PWD
parentdir="$(dirname "$dir")"
bindir=/bin
fullbindir=$parentdir$bindir

# Make the binary directory if it doesn't exist
if [ -d $fullbindir ]
  then
    echo "LSDTopoTools binary directory exists!"
  else
    echo "LSDTopoTools binary directory doesnt' exist. I'm making one"
    mkdir $fullbindir
fi


echo The full directory for binary files is:
echo $fullbindir

libdir=/lib
fulllibdir=$parentdir$libdir

# Make the binary directory if it doesn't exist
if [ -d $fulllibdir ]
  then
    echo "LSDTopoTools library directory exists!"
  else
    echo "LSDTopoTools library directory doesnt' exist. I'm making one"
    mkdir $fulllibdir
fi



cmake -DCMAKE_INSTALL_PREFIX=$PREFIX $SRC_DIR -DCMAKE_INSTALL_BINDIR=$fullbindir -DCMAKE_INSTALL_LIBDIR=$fullbindir

make
make install