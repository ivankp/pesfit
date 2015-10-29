CC  := gcc
CXX := g++

CXXFLAGS := -std=c++11 -Wall -O3

LIBS := -lboost_program_options

# Use boost libraries if GCC version doesn't support std equivalents
ifneq ($(shell printf "%s\n4.9\n" `gcc --version | sed -n 's/^gcc .* //p'` | sort -V | head -1),4.9)
	LIBS += -lboost_regex
endif

ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LIBS     := $(shell root-config --libs) -lMinuit2 -lRooFitCore -lRooFit -lRooStats

.PHONY: all clean

NODEPS := clean

SRCDIR := src
BLDDIR := .build
EXEDIR := bin

SRCS := $(shell find $(SRCDIR) -name "*.cc")
OBJS := $(patsubst $(SRCDIR)/%.cc,$(BLDDIR)/%.o,$(SRCS))
DEPS := $(patsubst %.o,%.d,$(OBJS))

GREP_EXE := grep -r '^ *int *main *(' $(SRCDIR) | cut -d':' -f1
EXES := $(patsubst $(SRCDIR)/%.cc,$(EXEDIR)/%,$(shell $(GREP_EXE)))

all: $(EXES)

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPS)
endif

# create dependencies
$(BLDDIR)/%.d: $(SRCDIR)/%.cc
	@echo DEP $(notdir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT '$(patsubst $(SRCDIR)/%.cc,$(BLDDIR)/%.o,$<)' $< -MF $@

# compile objects
$(BLDDIR)/%.o : $(SRCDIR)/%.cc
	@echo CXX $(notdir $@)
	@$(CXX) -c -I$(SRCDIR) $(CXXFLAGS) $(ROOT_CXXFLAGS) $< -o $@

# link executables
$(EXEDIR)/% : $(BLDDIR)/%.o
	@echo LD $(notdir $@)
	@$(CXX) $(filter %.o,$^) -o $@ $(LIBS) $(ROOT_LIBS)

# directories as order-only-prerequisites
$(OBJS) $(DEPS): | $(BLDDIR)
$(EXES): | $(EXEDIR)

# make directories
$(BLDDIR) $(EXEDIR):
	mkdir $@

clean:
	@rm -vfr $(BLDDIR) $(EXEDIR)

$(EXEDIR)/pesfit: $(BLDDIR)/TGraph_fcns.o
