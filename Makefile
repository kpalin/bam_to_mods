#
# Compiler flags
#
CC     := $(GXX)
CFLAGS += -Wall

#
# Project files
#
SRCS = bam_mods_to_text.cpp
OBJS = $(SRCS:.cpp=.o)
EXE  = bam_to_mods
LIBS :=   -lhts -lbz2 -lhts -lz -lm   -lpthread
#
# Debug build settings
#
DBGDIR = debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))  
DBGCFLAGS = $(DEBUG_CXXFLAGS) -g -O0 -DDEBUG  -fno-tree-vectorize  -Wextra # -fsanitize=address
DBGLDFLAGS = $(LDFLAGS)  
#
# Release build settings
#
RELDIR = release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = $(CXXFLAGS) -O3 -DNDEBUG

.PHONY: all clean debug prep release remake test testall

# Testing
BATS=./test/bats/bin/bats
TESTS=test/basics.bats

# Default build
all: prep release debug

test: $(DBGEXE) $(TESTS)
	$(BATS) $(TESTS)

testall: $(DBGEXE)
	$(BATS) test/

#
# Debug rules
#
debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CXX) $^  $(LDFLAGS) $(LIBS) $(DBGLDFLAGS)  -o $(DBGEXE)

$(DBGDIR)/%.o: src/%.cpp
	$(CXX) -c $(CFLAGS) $(DBGCFLAGS)  -o $@ $<

#
# Release rules
#
release: $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CXX)  $^ $(LDFLAGS) $(LIBS) $(LDFLAGS)  -o $(RELEXE)

$(RELDIR)/%.o: src/%.cpp
	$(CXX) -c $(CFLAGS) $(RELCFLAGS) -o $@ $<

#
# Other rules
#
prep:
	@mkdir -p $(DBGDIR) $(RELDIR)

remake: clean all

clean:
	rm -f $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS)
