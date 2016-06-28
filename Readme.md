# Readme.md  MJD_VetoTracks
#
#  David J Tedeschi - Univ. of South Carolina - 2015
#

Directories:

VetoDisplay - muon track visualization program
Landscape - Black Hills topology 
Data - all input dat files 
Doc -  Latex Documentation of this projecy


Building the APPs:
 
Descend into the VetoDisplay or Landscape directory and use the commands 
"make clean" then "make create" to compile and create an executable. 

Dependencies:
Must first have Root, MGDO, and GAT installed. 

Execution:
./Vetodisplay - run with no arguments reads all input files and generates pdf outputs

./VetoDisplay <runNumber> <eventCount>
If there is a certain muon event you'd like to see in a rotatable 3-D model, simply run from the VetoDisplay directory 
with the arguments above.


Output:
When run from the command line, VetoDisplay generates pdf files to the "output" folder. 
The subdirectories coorespond to number of panels hit, with 9+ panels hit going into the 'other' directory, 
as well as a final heatmap.

To clear all output subdirectories:
cd to the main output directory "output" and use the command:

find . ! -type d -exec rm '{}' \;

careful, dangourous command if used from any other directory!



Documentation:

For more information, see the Doc folder.

