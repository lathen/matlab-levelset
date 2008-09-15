
# Assume that Matlab binaries are in the path
MEX = mex
MEXEXT = mexext
CXXFLAGS= -I. -I../../MEX -O2 -fPIC

EXT = $(shell $(MEXEXT))
BIN_PATH = ..

TOOLS = private/diff_central_order2 \
        private/diff_upwind_order1 \
    	private/diff2_order2 \
        private/reinitialize_fastmarching \
        private/triangulate \
	private/min_curvature

TOOLS_MEX := $(addprefix $(BIN_PATH)/, $(TOOLS:=.$(EXT)))

all : $(TOOLS_MEX)

# Define how gcc makes mex-files
$(BIN_PATH)/%.$(EXT): %.cpp
	$(MEX) CXXFLAGS="$(CXXFLAGS)" $< -output $(basename $@)

# Automatic dependency updating
%.d: %.cpp
	@echo "Building dependencies for"  $(basename $@).$(EXT)
	@$(MEX) CXXFLAGS="$(CXXFLAGS) -MM -MT $(BIN_PATH)/$(basename $@).$(EXT)" -c $< > $@

ifneq ($(MAKECMDGOALS),clean)
include $(TOOLS:=.d)
endif

clean:
	rm -rf $(TOOLS_MEX)
	rm -rf $(TOOLS:=.d)
