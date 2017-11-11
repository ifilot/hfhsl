# set compiler and compile options
EXEC = hfcxx
CXX = g++                             # use the GNU C++ compiler
OPTS = -std=c++0x -O3 -Wall           # use some optimization, report all warnings and enable debugging
CFLAGS = $(OPTS)                      # add compile flags
LDFLAGS = # specify link flags here

# set a list of directories
OBJDIR = ./obj
BINDIR = ./bin
SRCDIR = ./src

# set the include folder where the .h files reside
CFLAGS += -I$(SRCDIR)

# specify OS-related stuff
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CFLAGS  += -I/opt/local/include/eigen3 -I/opt/local/include
	LDFLAGS +=
else
    CFLAGS += `pkg-config --cflags eigen3` -I/mingw64/include
	LDFLAGS +=
endif

# add here the source files for the compilation
SOURCES = hfcxx.cpp

# create the obj variable by substituting the extension of the sources
# and adding a path
_OBJ = $(SOURCES:.cpp=.o)
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

all: $(BINDIR)/$(EXEC)

$(BINDIR)/$(EXEC): $(OBJ)
	$(CXX) -o $(BINDIR)/$(EXEC) $(OBJ) $(LIBDIR) $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c -o $@ $< $(CFLAGS)

clean:
	rm -vf $(BINDIR)/$(EXEC) $(OBJ)
