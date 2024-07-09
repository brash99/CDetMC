name := AnaBarMC
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

cleanup:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_* *~ src/*~ include/*~ macros/*~
	rm -f *.so *.d

########################### ROOT #################################

ifdef ROOTSYS
ifndef G4UI_USE_ROOT
  ROOTCPPFLAGS   = $(shell $(ROOTSYS)/bin/root-config --cflags)
  CPPFLAGS      += -DG4ANALYSIS_USE_ROOT $(ROOTCPPFLAGS)
  ROOTLIBS       = $(shell $(ROOTSYS)/bin/root-config --nonew --glibs)
  ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
  ROOTLIBS      += -lMinuit
  LDLIBS        += $(ROOTLIBS)
endif
endif
