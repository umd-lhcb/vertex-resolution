BINPATH	:=	bin
VPATH	:=	src
CPP_FILES	:=	$(wildcard src/*.cpp)
EXE_FILES	:=	$(patsubst src/%.cpp,$(BINPATH)/%,$(CPP_FILES))

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags) -Iinclude
LINKFLAGS	:=	$(shell root-config --libs)
ADDCXXFLAGS	:=	-O2 -march=native -mtune=native


###########
# General #
###########

.PHONY: exe
exe: $(EXE_FILES)

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

$(BINPATH)/%: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS)
