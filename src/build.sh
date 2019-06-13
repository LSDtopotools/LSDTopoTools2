export LDFLAGS="-Wl,-rpath,$PREFIX/lib -L$PREFIX/lib $LDFLAGS"

# This sorts out the directory into which the code will be built
dir=$PWD
parentdir="$(dirname "$dir")"
bindir=/bin
fullbindir=$parentdir$bindir
echo The full directory for binary files is:
echo $fullbindir

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX $SRC_DIR -DCMAKE_INSTALL_BINDIR=$fullbindir

make
make install