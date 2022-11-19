BINPATH	:=	bin
VPATH	:=	src:$(BINPATH)

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags) -Iinclude
LINKFLAGS	:=	$(shell root-config --libs)


###########
# General #
###########

exe: PrintMCDecay ReweightRDX ReweightRDXDebug

.PHONY: clean
clean:
	@rm -rf ./bin/*
	@rm -rf ./gen/*


#########
# Plots #
#########


####################
# Generic patterns #
####################

.SECONDARY:

%: %.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS)
