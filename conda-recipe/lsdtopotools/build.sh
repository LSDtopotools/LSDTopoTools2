export LDFLAGS="-Wl,-rpath,$PREFIX/lib -L$PREFIX/lib $LDFLAGS"

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX $SRC_DIR -DCMAKE_INSTALL_BINDIR=bin

make
make install
