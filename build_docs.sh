#!/bin/bash
rm -rf html/
doxygen Doxyfile
cp img/LSD-logo.png html/
echo $?
