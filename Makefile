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
LIBS := -lhts $(LIBS) -lpthread
#
# Debug build settings
#
DBGDIR = debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))  
DBGCFLAGS = $(DEBUG_CXXFLAGS) -g -O0 -DDEBUG  -fno-tree-vectorize  -Wextra # -fsanitize=address
DBGLDFLAGS = $(LDFLAGS)   #-fsanitize=address 
#
# Release build settings
#
RELDIR = release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = $(CXXFLAGS) -O3 -DNDEBUG

.PHONY: all clean debug prep release remake

# Default build
all: prep release debug

#
# Debug rules
#
debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CXX) $(LIBS) $(DBGLDFLAGS)  -o $(DBGEXE) $^

$(DBGDIR)/%.o: src/%.cpp
	$(CXX) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $<

#
# Release rules
#
release: $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CXX) $(LIBS) $(LDFLAGS)  -o $(RELEXE) $^

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