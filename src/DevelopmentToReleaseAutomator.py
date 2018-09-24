
#===============================================================================
# This script is for automating the posting of the latest versions of our
# development code. It is a bridge between the private development library
# and our public github releases
# It 
# 1) Finds all the makefiles in a directory supplied by the user
# 2) Takes all the files that are used in the makefiles and copies them over
#    to the directory in the github folder
# 
#
# INSTRUCTIONS
# 1) You need the development code and the release code to have the same root folder. 
#    for example, I have development code in
#    /smudd/Git_projects/LSDTT_Development
#    and the release codes in
#    /smudd/Git_projects/NameOfRelease
# 2) Pick the code you want to release. This will be in a "driver_functions" folder.
# 3) Figure out the directory you want to put the release. 
# 4) In this code, search for HIBERNIANFC SCOTTISH CUP 2016 and jump there. 
#    Below that line you will find the direcories you need. 
# 5) Do not change ObjectsDirectory
# 6) DriverDirectory input the driver directory IN THE DEVELPMENT BRANCH
# 7) TargetDirectory input where the release repo is. It will be preceded by `../`
# 8) Run the python code. 
# 
#
#===============================================================================

import numpy as np
from glob import glob
import os
import shutil

# This function takes a string and looks for the seperators
# it then reformats to the correct operating system
def ReformatSeperators(path):
    
    path = RemoveEscapeCharacters(path)
    #print "Path is: " + path    
    
    # loook for the various seperators in the string
    if "\\" in path:
        #print "I found \\"
        path = RemoveEscapeCharacters(path)
        splitpath = path.split("\\")
    elif "/" in path:
        #print "I found ''/''"
        if "//" in path:
            #print "I found ''//''"
            splitpath = path.split("//")
            
        else:
            splitpath = path.split("/")
    else:
        print "I did not find a valid seperator"
        splitpath = "NULL"

    # Now reconstitute the path using the operating system seperator
    # get the size of the split path
    n_spaces = len(splitpath)   
    newpath = ""   
    for s in range (0,n_spaces-1):
        newpath = newpath+splitpath[s]+os.sep 
            
    # now the last element
    newpath = newpath+splitpath[n_spaces-1]
    
    newpath2 = RemoveEscapeCharacters(newpath)
    
    return newpath2    

# This function takes a path with any seperators and converts the seperators
# to the current operating system and adds a seperator at the end of the
# path                    
def AppendSepToDirectoryPath(path):
    
    #first check to see if there are escape characters
    RemoveEscapeCharacters(path)    
    
    # now reformat the seperators
    newpath = ReformatSeperators(path)
    
    # now check to see if the last character is a sep
    if newpath[-1] != os.sep:
        newpath = newpath+os.sep
    return newpath

# This function takes a filename (that could include a full directory path)
# and lops off the path as well as the extension
def GetFilePrefix(filename):
    newfilename = ReformatSeperators(filename)
    
    # split the file
    splitname = newfilename.split(os.sep)
    fname = splitname[-1]
    splitfname = fname.split(".")
    
    fileprefix = splitfname[0]
    return fileprefix

# This gets the filename without the path
def GetFileNameNoPath(filename):
    newfilename = ReformatSeperators(filename)
    
    # split the file
    splitname = newfilename.split(os.sep)
    fname = splitname[-1]
    return fname


# This gets the last directory level
# say if the path is home\yo\ma\
# then this function returns ma
def GetLastDirectoryLevel(path):
    newpathname = ReformatSeperators(path)
    pathname = AppendSepToDirectoryPath(newpathname)
    
    # now split the path
    splitpath = pathname.split(os.sep)
    
    # now remove the final path level
    n_spaces = len(splitpath)   
    newpath = splitpath[n_spaces-2] 

    return newpath
        
        
# This gets the path of a filename (with a full directory path)
def GetPath(filename):
    newfilename = ReformatSeperators(filename)
    
    splitname = newfilename.split(os.sep)
    n_levels = len(splitname)   
    newpath = ""   
    for s in range (0,n_levels-1):
        newpath = newpath+splitname[s]+os.sep
        
    return newpath

# This function finds the level of the path
def GetPathLevel(path):
    newpathname = ReformatSeperators(path)
    pathname = AppendSepToDirectoryPath(newpathname)
    splitname = pathname.split(os.sep)
    n_levels = len(splitname)  
    return n_levels
    
    
# This removes one level from the directory    
def RemoveDirectoryLevel(path):
    # Make sure the formatting is appropriate and there is a seperator at the end
    newpathname = ReformatSeperators(path)
    pathname = AppendSepToDirectoryPath(newpathname)
    
    # now split the path
    splitpath = pathname.split(os.sep)
    
    # now remove the final path level
    n_spaces = len(splitpath)   
    newpath = "" 
    for s in range (0,n_spaces-2):
        newpath = newpath+splitpath[s]+os.sep
        
    return newpath
    


# This is necessary because of stupid windows seperators    
def RemoveEscapeCharacters(line):
    line = line.rstrip().replace('\n', '\\n')
    line = line.rstrip().replace('\b', '\\b') 
    line = line.rstrip().replace('\r', '\\r')  
    line = line.rstrip().replace('\t', '\\t')  
    line = line.rstrip().replace('\n', '\\n')  
    line = line.rstrip().replace('\f', '\\f')  
    line = line.rstrip().replace('\v', '\\v') 
        
    # this last one deals with the infuriating special case of \b
    line = line.rstrip().replace('\x08', '\\b')     
    
    # get rid of leading and trailing spaces
    line = line.strip()     
    
    return line
    
def RemoveWhitespace(line):
    line = line.replace(" ", "")
    return line



# This goes into a data directory and trawls any makefiles for the relevant .cpp and .hpp
# files necessary to compile the code. 
def GetRequiredFilesFromFolder(DataDirectory):
    
    #print "Current directory is: " + os.getcwd()   
    
    #Append a closing slash to the data directory if one is not already there
    NewDataDirectory = ReformatSeperators(DataDirectory)   
    DataDirectory = AppendSepToDirectoryPath(NewDataDirectory)
    
    #print "DataDirectory is (2): " + DataDirectory    

    # Start with a master list of required files
    required_files = []

    # find all the makefiles in the directory
    for FileName in glob(DataDirectory+"*.make"):
        #print "FileName is: " + FileName
        
        # Now you need to find all the sources in the makefile
        f = open(FileName,'r')  # open file
        lines = f.readlines()   # read in the data
        f.close()
        
        # Initiate an empty list that will contain the filenames
        cppfiles = []
        
        # add the makefile to the files to be copied
        required_files.append(FileName)
        
        # loop through the lines, flag a start where SOURCES line starts
        start_flag = 0
        for line in lines:
            if "SOURCE" in line:
                start_flag = 1
            
            # If there is OBJECTS in the line, stop looking for ".cpp"
            if "OBJECTS" in line:
                start_flag = 0
            
            # Look for .cpp between SOURCES and OBJECTS
            if start_flag == 1:
                if ".cpp" in line:
                    # seperate the line using spaces
                    split_line = line.split(' ')
                    this_item = ""
                    
                    # Go through the split line looking for .cpp
                    for item in split_line:
                        if ".cpp" in item:
                            this_item = item
                            
                            # get rid of the SOURCES
                            new_this_item = this_item.replace("SOURCES=","")
                            
                    #print "This cpp file is: " + new_this_item
                    
                    # get rid of stupid escape characters
                    this_file = RemoveEscapeCharacters(new_this_item)
                    
                    cppfiles.append(this_file)
                    
        # now print to screen the files required for this makefile
        #print "The files required for this makefile are: "
        #print cppfiles
        
        # now append to directory...this requires some logic because of the ../ seperators
        for filename in cppfiles:
            
            #print "Filename is: " + filename            
            
            # special logic for the ../ seperator            
            if "../" in filename:
                #print "There is a lower level in this filename, this means it is an object"
                thisfile =  filename.replace("../","")
                thisdirectory = RemoveDirectoryLevel(DataDirectory)
                fullfile = thisdirectory+thisfile
                
                fullfile2 = fullfile.replace(".cpp",".hpp")
                required_files.append(fullfile2)
            else:
                fullfile = DataDirectory+filename  
        
            # append to the required files list
            required_files.append(fullfile)
                
    # now thin out the required files to remove duplicates
    nd = set(required_files)
    required_files_noduplicates = list(nd)

    #print "/n/n================================="    
    #print "Required files are: " 
    #print required_files
    #print "--------"
    #print "And removing duplicates:"
    #print required_files_noduplicates
    #print "====================================="
    
    return required_files_noduplicates

# This function checks the file structures and either makes directories or
# throws errors when file structures do not exist
def CheckFileStructuresForCopy(ObjectsDirectory,DriverDirectory,TargetDirectory):
    # Format the target directories
    Td = ReformatSeperators(TargetDirectory)   
    TargetDirectory = AppendSepToDirectoryPath(Td)  
    TDd = TargetDirectory + DriverDirectory
    TargetDriverDirectory = AppendSepToDirectoryPath(TDd)

    # Format the source directories
    Od = ReformatSeperators(ObjectsDirectory)   
    ObjectsDirectory = AppendSepToDirectoryPath(Od)
    Dd = ObjectsDirectory+DriverDirectory
    DriverDirectory = AppendSepToDirectoryPath(Dd)
            
    # Check if the source directories exist
    if not os.access(ObjectsDirectory,os.F_OK):
        print "The object directory for the code doesn't exist!"
        print "You wanted this directory: " + ObjectsDirectory
        return 0
    if not os.access(ObjectsDirectory,os.F_OK):
        print "The driver directory for the code doesn't exist!"
        print "You wanted this directory: " + DriverDirectory
        return 0        
    if not os.access(ObjectsDirectory+"TNT"+os.sep,os.F_OK):
        print "The TNT directory for the code doesn't exist!"
        print "You wanted this directory: " + ObjectsDirectory+"TNT"+os.sep
        return 0 
     
    # check if the target object directory exists
    if not os.access(TargetDirectory,os.F_OK):
        print "The target directory for the code doesn't exist!"
        print "You wanted this directory: " + TargetDirectory
        print "I am making that now, along with the driver directory"
        os.mkdir(TargetDirectory)
        if not os.access(TargetDirectory,os.F_OK):
            print "WTF the directory was not made??!"
        os.mkdir(TargetDriverDirectory)
        
    # check just the driver directory
    if not os.access(TargetDriverDirectory,os.F_OK):
        print "The target driver directory for the code doesn't exist!"
        print "You wanted this directory: " + TargetDriverDirectory
        print "I am making that now"
        os.mkdir(TargetDriverDirectory)    
        
    # Check if the TNT directory exists. If it does, remove and replace it
    # If it doesn't , just copy it across
    TNTTargetDirectory = TargetDirectory+'TNT'+os.sep
    TNTSourceDirectory = ObjectsDirectory+'TNT'+os.sep
    if not os.access(TNTTargetDirectory,os.F_OK):
        print "The target TNT directory for the code doesn't exist!"
        print "You wanted this directory: " + TargetDriverDirectory
        print "I am making that now"
        shutil.copytree(TNTSourceDirectory,TNTTargetDirectory)
    else:
        print "There is a TNT directory here already. Removing and replacing"
        shutil.rmtree(TNTTargetDirectory)
        shutil.copytree(TNTSourceDirectory,TNTTargetDirectory)

    print "========================="
    print "DriverDirectory: " + DriverDirectory 
    print "ObjectsDirectory: " + ObjectsDirectory
    print "TargetDirectory: " + TargetDirectory
    print "TargetDriverDirectory: " + TargetDriverDirectory
    print "========================="
        
    return ObjectsDirectory,DriverDirectory,TargetDirectory,TargetDriverDirectory 
    

# This function is for copying a group of files from the makefile in a driver directory
def CopyRequiredFilesToGitRepository(ObjectsDirectory,DriverDirectory,TargetDirectory):

    # Ensure the directories exist
    ObjectsDirectory,DriverDirectory,TargetDirectory,TargetDriverDirectory = CheckFileStructuresForCopy(ObjectsDirectory,DriverDirectory,TargetDirectory)
             
    # Now get the required files
    print "\n\n\n================================="
    required_files_noduplicates = GetRequiredFilesFromFolder(DriverDirectory) 
    print "The required files are: "
    print required_files_noduplicates
    print "================================="


    # loop through these files, collecting the filenames and directory names
    # first you need to know what directory level the driver files are in
    print "\n\n\n======================================"
    n_level_of_driver_directory = GetPathLevel(DriverDirectory)
    for FileName in required_files_noduplicates:
        # you need to know what level the file is
        ThisPath = GetPath(FileName)
        ThisLevel = GetPathLevel(ThisPath)
        
        #if it is the same level as the driver directory, it is in the driver directory!
        if ThisLevel == n_level_of_driver_directory:        
            CopyDirectory = TargetDriverDirectory
            CopyFileName = GetFileNameNoPath(FileName)
            CopyFileNameWithPath = CopyDirectory+CopyFileName
        else:
            CopyDirectory = TargetDirectory
            CopyFileName = GetFileNameNoPath(FileName)
            CopyFileNameWithPath = CopyDirectory+CopyFileName            
            
        print "The filename is: " + FileName
        print "The copy filename is: " + CopyFileNameWithPath
        shutil.copy(FileName, CopyFileNameWithPath)
        
        # now copy the files over
        
        
    print "=============================================="             
                    
# This is to test if the a make file actually results in a program
def TestMake(Path,Makefile):
    """
    Make the file Makefile stored in Path and return True if the returncode is
    0 (Success). If the returncode is 2 (Failure) or any other value (wierd behaviour) 
    return False.
    
    Gives no feedback on why the make failed. Not tested on Windows.
    
    SWDG 26/8/15
    """
    
    from subprocess import Popen, PIPE    
    
    cmd = 'cd '+Path + '; make -f '+Makefile
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)    
    #subprocess.call(cmd, shell=True)
    p.communicate()[0]
    returncode = p.returncode 
    
    if returncode == 0:
        return True
    else:
        return False
    

if __name__ == "__main__":
    # YOU NEED TO MODIFY THIS DIRECTORY
    
    # This one is for running in windows
    # If you are in windows uncomment and comment out the linux version
    #ObjectsDirectory = 'T:\devel_projects\LSDTopoTools\trunk'
    #DriverDirectory = 'Analysis_driver'
    #TargetDirectory = 'T:\Git_projects\LSDTopoTools_CRNBasinwide'
    
    # This one is for running directly in linux
    # If you are in linux uncomment and comment out the windows version
    # Search for HIBERNIANFC SCOTTISH CUP 2016 to get to the correct line 
    ObjectsDirectory = './'
    DriverDirectory = './driver_functions_MuddChi2014'
    TargetDirectory = '../LSDTopoTools_ChiMudd2014/'   
    #DriverDirectory = './driver_functions_ChannelExtraction'
    #TargetDirectory = '../LSDTopoTools_ChannelExtraction/'    
    #DriverDirectory = './Analysis_driver'
    #TargetDirectory = '../LSDTopoTools_AnalysisDriver/'  
    
    CopyRequiredFilesToGitRepository(ObjectsDirectory,DriverDirectory,TargetDirectory)    
    
    
    #DataDirectory =  'T:\devel_projects\LSDTopoTools\trunk\driver_functions_MuddChi2014'
    # THis one is for running directly in linux    
    #DataDirectory =  '/home/smudd/SMMDataStore/devel_projects/LSDTopoTools/trunk/driver_functions_MuddChi2014'        
    #required_files_noduplicates = GetRequiredFilesFromFolder(DataDirectory)   
    
    #print required_files_noduplicates
