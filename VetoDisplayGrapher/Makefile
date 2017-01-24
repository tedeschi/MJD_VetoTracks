include $(MGDODIR)/buildTools/config.mk
include $(MGDODIR)/buildTools/BasicROOTMakefile

# Give the application name, which must be the stem of cc file with 'main'.  
APPS = VetoDisplayGrapher
# The next two lines are important
SHLIB = 
ARCHIVE = 
TAMDIR ?= $(ROOTSYS)
# Include the correct flags, 
INCLUDEFLAGS = $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms -I$(MGDODIR)/Tabree -I$(MGDODIR)/MJDB
INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(ROOTSYS)/include -I$(TAMDIR)/inc -I$(TAMDIR)/include
INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis
LIBFLAGS = -L$(ROOTSYS)/lib -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOTabree -lMGDOMJDB -lGeom -lEG -lMathMore
	
all: create
	
create: VetoDisplayGrapher.o
	g++ -o VetoDisplayGrapher VetoDisplayGrapher.o $(INCLUDEFLAGS) $(LIBFLAGS) `root-config --libs`

VetoDisplayGrapher.o: VetoDisplayGrapher.cc VetoDisplayGrapher.hh
	g++ -o VetoDisplayGrapher.o VetoDisplayGrapher.cc -c $(INCLUDEFLAGS) $(LIBFLAGS) `root-config --cflags`
