BINPATH	:=	bin
VPATH	:=	src
CPP_FILES	:=	$(wildcard src/*.cpp)
EXE_FILES	:=	$(patsubst src/%.cpp,$(BINPATH)/%,$(CPP_FILES))

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags) -Iinclude
LINKFLAGS	:=	$(shell root-config --libs)
ADDCXXFLAGS	:=	-O2 -march=native -mtune=native

OS := $(shell uname)
ifeq ($(OS),Darwin)
  $(info OS is $(OS) (macOS), adding -lc++fs to LINKFLAGS)
  LINKFLAGS := $(shell root-config --libs) -lc++fs
endif

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


########
# Test #
########
.PHONY: test-all

test-all: gen/test-rdx-run2-Bd2Dst0MuNu-sim09k.root

gen/test-%.root: samples/%.root $(BINPATH)/ApplyVertexSmear
	@./bin/ApplyVertexSmear -i $< -x ./inputs/smearing_vec.root -o $@


####################
# Generic patterns #
####################

$(BINPATH)/%: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS)
