# MJVetoTracks
Decend into the VetoDisplay or Landscape directory and use the commands "make clean" then "make create" to compile and create an executable. Must first have Root and MGDO installed. When run from the command line, VetoDisplay generates pdf files to the "output" folder. The subdirectories coorespond to number of panels hit, with 9+ panels hit going into the 'other' directory, as well as a final heatmap.

If there is a certain run you'd like to see a rotatable 3-D model of, simply run from the VetoDisplay directory: VetoDisplay runNumber eventCount

For more information, see the Doc folder.

To clear all output subdirectories:
cd to the main output directory "output" and use the command:

find . ! -type d -exec rm '{}' \;

careful, dangourous command if used from any other directory!

