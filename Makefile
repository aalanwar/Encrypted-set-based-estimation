BINDBG := HW_debug
BIN := HW
MAINSRC := amr_test.cpp 
OBJDIR := .o
OBJDBGDIR := .od
DEPDIR := .d
CXX := g++
CXXFLAGS := --std=gnu++17  -I/usr/include/eigen3 -I/usr/local/MATLAB/R2019b/extern/include -Irepresentations -Icrypto   
LIBS := -L/usr/local/MATLAB/R2019b/bin/glnxa64 -Lcrypto -Lcrypto/math
SRCS     :=                \
	representations/eigen2mat/eigen2mat.cpp \
	representations/CSVReader.cpp \
	representations/EncZonotope/EncZonotope.cpp \
	representations/EncConZonotope/EncConZonotope.cpp \
	representations/EncStrip/EncStrip.cpp \
	representations/Strip/Strip.cpp \
	representations/EncEntity/EncEntity.cpp \
	representations/Zonotope/Zonotope.cpp \
	representations/ConZonotope/ConZonotope.cpp \
	$(wildcard crypto/*.cc) \
	$(wildcard crypto/math/*.cc)
OBJS := $(patsubst %,$(OBJDIR)/%.o,$(basename $(SRCS)))
OBJSDBG := $(patsubst %,$(OBJDBGDIR)/%.o,$(basename $(SRCS)))
DEPS := $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

$(shell mkdir -p $(dir $(OBJS)) >/dev/null)
$(shell mkdir -p $(dir $(OBJSDBG)) >/dev/null)
$(shell mkdir -p $(dir $(DEPS)) >/dev/null)


LINK := $(CXX) $(MAINSRC) $(CXXFLAGS) -Wl,-rpath=/usr/local/MATLAB/R2019b/bin/glnxa64 $(OBJS) -o $@ 

LINKDBG := $(CXX) -g $(MAINSRC) $(CXXFLAGS) -Wl,-rpath=/usr/local/MATLAB/R2019b/bin/glnxa64 $(OBJSDBG) -o $@ 

COMPILE:= $(CXX) $(CXXFLAGS) $(LIBS) -lntl -lgmp -lgmpxx -lpthread -lmat -lmx  -c -o $@

COMPILEDBG:= $(CXX) -g $(CXXFLAGS) $(LIBS) -lntl -lgmp -lgmpxx -lpthread -lmat -lmx -c -o $@

all: $(BIN) 

debug: $(BINDBG)

clean: 
	rm -f -r $(OBJDIR) $(OBJDBGDIR) $(DEPDIR) 

# BIN is rebuilt when any OBJS changes or MAINSRC
# HOW to build BIN --> using the second line that includes LINK
$(BIN): $(OBJS) $(MAINSRC)
		$(LINK) $@ $(LIBS) -lntl -lgmp -lgmpxx -lpthread -lmat -lmx 

#Same as before
$(BINDBG): $(OBJSDBG) $(MAINSRC)
		   $(LINKDBG) $@ $(LIBS) -lntl -lgmp -lgmpxx -lpthread -lmat -lmx 

#for each cpp file generates .o file 
$(OBJDIR)/%.o: %.cpp
			$(COMPILE) $@ $<

#for each cc file generates .o file 
$(OBJDIR)/%.o: %.cc
			$(COMPILE) $@ $<

$(OBJDBGDIR)/%.o: %.cpp
			$(COMPILEDBG) $@ $<

$(OBJDBGDIR)/%.o: %.cc
			$(COMPILEDBG) $@ $<					

.PRECIOUS: $(DEPDIR)/%.d 
$(DEPDIR)/%.d: ;

-include $(DEPS)

