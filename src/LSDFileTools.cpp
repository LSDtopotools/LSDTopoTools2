
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "LSDFileTools.hpp"
using namespace std;


#ifndef FileTools_CPP
#define FileTools_CPP

bool DoesFileExist(string filename)
{
  bool WellDoesIt = false;
  
  struct stat sb;

  if (stat(filename.c_str(), &sb) == 0 && S_ISREG(sb.st_mode))
  {
    WellDoesIt = true;
  }
  return WellDoesIt;
}

bool DoesFileExist(string path, string fname)
{
  bool WellDoesIt = false;
  
  string filename = path+fname;
  
  struct stat sb;

  if (stat(filename.c_str(), &sb) == 0 && S_ISREG(sb.st_mode))
  {
    WellDoesIt = true;
  }
  return WellDoesIt;
}

bool DoesDirectoryExist(string dirname)
{
  bool WellDoesIt = false;

  struct stat sb;

  if (stat(dirname.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
  {
    WellDoesIt = true;
  }
  return WellDoesIt;
}


#endif
