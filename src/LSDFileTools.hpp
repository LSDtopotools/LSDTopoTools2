
#include <string>
using namespace std;


#ifndef FileTools_HPP
#define FileTools_HPP

// These functions check to see if a file or a directory exists
bool DoesFileExist(string filename);
bool DoesFileExist(string path, string fname);
bool DoesDirectoryExist(string dirname);


#endif
