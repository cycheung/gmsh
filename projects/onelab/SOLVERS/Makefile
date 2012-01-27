SOLVERS=onelab
SOURCES=OS StringUtils OnelabClients OnelabMessage main

COMMON_DIR=$(GMSH_DIR)/Common
ONELAB_DIR=$(GMSH_DIR)/projects/onelab

OBJDIR=objects

ifeq ($(MINGW), 32)
	CXX=i686-w64-mingw32-g++
	LDFLAGS+=-static-libgcc -static-libstdc++ -static -lws2_32
	OBJDIR=objects-mingw32
	EXESUFFIX=32.exe
endif
ifeq ($(MINGW), 64)
	CXX=x86_64-w64-mingw32-g++
	LDFLAGS+=-static-libgcc -static-libstdc++ -static -lws2_32
	OBJDIR=objects-mingw64
	EXESUFFIX=64.exe
endif

TARGETS=$(SOLVERS:=_client$(EXESUFFIX))
all : $(TARGETS)


CPPFLAGS+=-I$(COMMON_DIR) -I$(ONELAB_DIR)
VPATH=$(ONELAB_DIR):$(COMMON_DIR):$(ONELAB_DIR)/SOLVERS

OBJECTS=$(foreach OBJ, $(SOURCES), $(OBJDIR)/$(OBJ).o)
$(OBJDIR)/%.o : %.cpp
	mkdir -p $(OBJDIR) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%_client$(EXESUFFIX) : $(OBJDIR)/%.o $(OBJECTS)
	$(CXX) $< $(OBJECTS) $(LDFLAGS) -o $@

clean :
	$(RM) $(TARGETS) -r $(OBJDIR)

.PRECIOUS: $(OBJDIR)/%.o
