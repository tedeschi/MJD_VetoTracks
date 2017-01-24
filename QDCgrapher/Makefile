include $(MGDODIR)/buildTools/config.mk
include $(MGDODIR)/buildTools/BasicROOTMakefile

# Give the application name, which must be the stem of cc file with 'main'.  
APPS = grapher
# The next two lines are important
SHLIB = 
ARCHIVE = 
TAMDIR ?= $(ROOTSYS)
# Include the correct flags, 
INCLUDEFLAGS = $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms -I$(MGDODIR)/Tabree -I$(MGDODIR)/MJDB
INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(ROOTSYS)/include -I$(TAMDIR)/inc -I$(TAMDIR)/include
INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis
LIBFLAGS = -L$(ROOTSYS)/lib -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOTabree -lMGDOMJDB -lGeom -lEG
	
all: create
	
#create: VetoDisplay.o
#	g++ -o VetoDisplay VetoDisplay.o $(INCLUDEFLAGS) $(LIBFLAGS) `root-config --libs`

#VetoDisplay.o: VetoDisplay.cc VetoDisplay.hh
#	g++ -o VetoDisplay.o VetoDisplay.cc -c $(INCLUDEFLAGS) $(LIBFLAGS) `root-config --cflags`

create: grapher.o
	g++ -o grapher grapher.o $(INCLUDEFLAGS) $(LIBFLAGS) `root-config --libs`

grapher.o: grapher.cc grapher.hh
	g++ -o grapher.o grapher.cc -c $(INCLUDEFLAGS) $(LIBFLAGS) `root-config --cflags`
