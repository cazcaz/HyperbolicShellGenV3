# HyperbolicShellGenV3
Not actually the third version, this is just the second model.

1) To install clone the repository to a file somewhere on a linux system (to use on windows, file paths in the code would need to be changed as they were not set up to be general)

2) install these modules with "sudo apt install " then the module name:
pip
python3
eigen
cmake

3) install these python modules with "sudo pip3 install " then the module name:
matplotlib
--upgrade matplotlib (Can't remember what this one did but I ran it)
numpy-stl
scipy

4) Install a c++ compiler, I used GCC 11.3.0 x86_64 

5) Build the executable with cmake (I did this through vscode using the cmake tools extension)

6) Run by running "sudo ./build/ShellGeneration" (maybe not with sudo if are concered, but with enough admin privileges to allow the program to edit and delete folders)

TO SETUP MODEL PARAMETERS

1) Go into the main.cpp file, and follow the comments
2) Create a shellParams struct, and set the desired parameters before pushing it to a list of shellParameters (To view which parameters are available simply open shellParams.h)
3) Strongly reccomend setting surface index to be unique for each surface so that no folder filenames are the same
4) Call a batchGen object to calculate them all, this will count how many compute threads are available and use every single one (if you dont want this then go into the batchGen class and change the m_threadCount after it has been initialised, i.e. take away 1 or 2 so it wont use every one)
5) Build and run executable, and outputs will be placed in separate folders in the surfaces folder in the repository

NOTE: When the code is run all previous surface folders are deleted and the new ones are created. The python code run at the end of main.cpp searches through every folder to find .txt files and converts them into graph .png and .stl model files placed in the respective surface folders

RECCOMENDATIONS ON PARAMETERS:
1) Dont change springcoeff (It likes being large)
2) Anything at 100 or below bending stiffness (circum) seems to fail surface quickly
3) Anything too small or big for the radius doesnt like existing
