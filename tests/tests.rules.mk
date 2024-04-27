# # Example usage:
# # 
# # 
# # include ../../common.compile.mk 
# # 
# # include ../../common.upkeep.mk
# # 
# # depdir := .deps
# # 
# # context=dense
# # contextdir=.
# # affix.$(context)=dense sparse operators combinatorics
# # projectdir=../../
# # pathvar=$(shell pwd)/../../


################################################################################
# EXPECTED VARIABLES:
# - projectdir : path/to/project/directory
# - context    : path/to/test/component/directory
# - contextdir : name of the test component
ifndef projectdir
$(error Expect 'projectdir')
endif
ifndef context
$(error Expect 'context')
endif
ifndef contextdir
$(error Expect 'contextdir')
endif

$(context).%: mycontext    := $(context)
$(context).%: mycontextdir := $(contextdir)

######################################################################################
# Decide whether to use .exe or .out as executable file ending 

ifeq ($(OS),Windows_NT)
ending := exe
else
ending := out
endif



######################################################################################
# Set the variables for this file 
# determine whether to use static or dynamic linking 

$(context).sources := $(sort $(wildcard $(contextdir)/*.cpp))

$(context).outs    := $(patsubst %.cpp,%.$(ending),$($(context).sources))

$(context).depdir  := $(contextdir)/.deps

$(context).dependencies := $(patsubst $(contextdir)/%.cpp,$(contextdir)/.deps/%.d,$($(context).sources))

$(context).include := $(patsubst %,-L$(projectdir)/%,$(affix.$(context)))
linkerprefix       :=-Wl,
$(context).rpath_t := $(patsubst %,-rpath=$(pathvar)/%,$(affix.$(context))) 
$(context).rpath   := $(patsubst %,$(linkerprefix)%,$($(context).rpath_t)) 

$(context).olib    := $(foreach X, $(affix.$(context)), $(projectdir)/$(X)/lib$(X).o  )
$(context).alib    := $(foreach X, $(affix.$(context)), $(projectdir)/$(X)/lib$(X).a  )
$(context).solib   := $(foreach X, $(affix.$(context)), $(projectdir)/$(X)/lib$(X).so )

# $(context).olib    := $(patsubst %,$(projectdir)/%/.all.o, $(affix.$(context)))
# $(context).alib    := $(patsubst %,$(projectdir)/%/lib%.a, $(affix.$(context)))
# $(context).solib   := $(patsubst %,$(projectdir)/%/lib%.so,$(affix.$(context)))

$(context).linklib     := $(patsubst %,-l%,$(affix.$(context)))
$(context).linkalib    := $(patsubst %,-l:lib%.a,$(affix.$(context)))
$(context).linksolib   := $(patsubst %,-l:lib%.so,$(affix.$(context)))

$(context).relevantlibs := 


ifeq ($(LINKINGTYPE),static)

$(context).mylib := $($(context).linkalib)
$(context).relevantlibs := $($(context).alib)

else ifeq ($(LINKINGTYPE),dynamic)
$(context).mylib := $($(context).linksolib)
#$(context).relevantlibs := 

else ifeq ($(LINKINGTYPE),objectfile)

$(context).mylib := $($(context).olib)
$(context).relevantlibs := $($(context).olib)

else ifeq ($(LINKINGTYPE),unspecified)

$(context).mylib := $($(context).linklib)
#$(context).relevantlibs := 

else

$(error No linking mode recognized: $(LINKINGTYPE)) 

endif


###################################################################################################
# All executables compiled depend on the makefiles 

%.$(ending): $(testsdir)/makefile $(testsdir)/tests.rules.mk $(testsdir)/tests.affices.mk $(projectsdir)/common.compile.mk

##########################################################################
# Specific definitions and recipes on how to compile the executables 

.PHONY: $($(context).depdir)
$($(context).depdir):
	@"mkdir" -p $@

DEPFLAGS = -MT $@ -MF $($(mycontext).depdir)/$*.d -MP -MMD
# We generate dependency files during compilation using the following compiler flags 
# -MT $@ : sets the target of the makefile rule, stripped of any directory components
# -MP    : add dependencies as phony targets
# -MF ...: sets the output file for the rules 
# -MMD   : list headers as a by-product of compiling, excluding system headers


$($(context).outs): mycontext    := $(context)
$($(context).outs): mycontextdir := $(contextdir)
$($(context).outs): $(contextdir)/%.$(ending): $(contextdir)/%.cpp $($(context).relevantlibs) | $($(context).depdir)
#	@ echo 
#	@ echo link type:   $(LINKINGTYPE)
#	@ echo context:     $(mycontext)
#	@ echo context dir: $(mycontextdir)
#	@ echo target:      $@
#	@ echo source file: $<
#	@ echo all prereq:  $^
#	@ echo contextdir:  $(mycontextdir)
#	@ echo depdir:      $($(mycontext).depdir)
#	@ echo include:     $($(mycontext).include)
#	@ echo rpath:       $($(mycontext).rpath)
#	@ echo object:      $($(mycontext).olib)
#	@ echo a:           $($(mycontext).alib)
#	@ echo so:          $($(mycontext).solib)
#	@ echo lib:         $($(mycontext).linklib)
#	@ echo a:           $($(mycontext).linkalib)
#	@ echo so:          $($(mycontext).linksolib)
#	@ echo used lib:    $($(mycontext).mylib)
	@echo Compiling $@ ...
#	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -std=c++2a -MT $@ -MF $($(mycontext).depdir)/$*.d -MM $(mycontextdir)/$*.cpp
ifeq ($(LINKINGTYPE),dynamic)
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include) $($(mycontext).rpath) $($(mycontext).mylib) -o $@ $(LDFLAGS) $(LDLIBS) $(DEPFLAGS)
else ifeq ($(LINKINGTYPE),static)
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include)                       $($(mycontext).mylib) -o $@ $(LDFLAGS) $(LDLIBS) $(DEPFLAGS)
else ifeq ($(LINKINGTYPE),objectfile)
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include)                       $($(mycontext).mylib) -o $@ $(LDFLAGS) $(LDLIBS) $(DEPFLAGS)
else 
	@echo Unable to compile the final object $@!
endif
	@touch $@
	
$($(context).dependencies):

-include $($(context).dependencies)

$($(context).outs): $(contextdir)/%.$(ending): $(contextdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/*.mk $(projectdir)/tests/*.mk

.PHONY: $(context).tests
$(context).tests: $($(context).outs)





###################################################
# How to automatically run the executables

$(context).runs        := $(patsubst %.cpp,%.run,$($(context).sources))

run: $(context).run
$(context).run: $($(context).runs)
$($(context).runs): %.run : %.$(ending)
	./$< 

$(context).silent_runs := $(patsubst %.cpp,%.silent_run,$($(context).sources))

silent_run: $(context).silent_run
$(context).silent_run: $($(context).silent_runs)
$($(context).silent_runs): %.silent_run : %.$(ending)
	./$< > /dev/null 

.PHONY: run silent_run $(context).run $(context).silent_run $($(context).runs) $($(context).silent_runs)

# # 2> /dev/null






# TODO: clean out the stuff below and adapt to the test directory


########################################################################
# Check whether the source files have correct syntax. Read-only.

$(context).sourcechecks := $(patsubst %.cpp,check-%.cpp,$($(context).sources))

.PHONY: checksources $($(context).sourcechecks)
checksources: $($(context).sourcechecks)
$($(context).sourcechecks): check-%.cpp : 
	$(info Check source: $*.cpp)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $*.cpp -fsyntax-only




########################################################################
# Apply clang-tidy to all cpp and hpp files in the directory. Read-only.

.PHONY: tidy $(context).tidy
tidy: $(context).tidy
$(context).tidy:
	clang-tidy $(mycontextdir)/*.?pp --config-file=$(projectdir)/.Tools/clang-tidy.yaml -- -std=c++17 -fno-exceptions


########################################################################
# Apply cppcheck to all cpp and hpp files in the directory. Read-only. 

.PHONY: cppcheck $(context).cppcheck
cppcheck: $(context).cppcheck
$(context).cppcheck:
	cppcheck -i ./.playground/ -i ./.legacy/ -i ./external/ \
	--enable=warning,style,performance,portability \
	--suppress=duplicateCondition --suppress=toomanyconfigs --suppress=sizeofFunctionCall --suppress=variableScope --suppress=missingReturn --suppress=assertWithSideEffect --suppress=useStlAlgorithm --suppress=knownConditionTrueFalse --suppress=unsignedPositive \
	--std=c++17 -q $(mycontextdir)/*pp



########################################################################
# Regex several useful things. Read-only. 
# - find trailing white spaces 
# - find non-ASCII characters 
# - find consecutive spaces 

.PHONY: grepissues $(context).grepissues
grepissues: $(context).grepissues
# $(context).grepissues: mycontext := $(context)
# $(context).grepissues: mycontextdir := $(contextdir)
$(context).grepissues:
#	@echo Search trailing whitespace...
#	@-grep --line-number --color '\s+$$' -r $(mycontextdir)/*pp
#	@echo Search non-ASCII characters...
#	@-grep --line-number --color '[^\x00-\x7F]' -r $(mycontextdir)/*pp
#	@echo Find consecutive spaces...
#	@-grep --line-number --color '\b\s{2,}' -r $(mycontextdir)/*pp
#	@echo Find standard asserts...
#	@-grep --line-number --color 'assert(' $(mycontextdir)/*pp
#	@echo Find usage of 'cout' ...
	@-grep --line-number --color 'cout' $(mycontextdir)/*pp
#	@echo Find floating-point numbers ...
#	@-grep --line-number --color -E '\.*[0-9]' $(mycontextdir)/*pp
#	@-grep --line-number --color -E '(0-9)e' $(mycontextdir)/*pp
	@-grep --line-number --color -E '([0-9]+e[0-9]+)|([0-9]+\.[0-9]+)|((+-\ )\.[0-9]+)|((+-\ )[0-9]+\.)' $(mycontextdir)/*pp


########################################################################
# Target 'check' is a generic test. Currently, it defaults to 'tidy'

check: tidy

.PHONY: check


########################################################################
# Commands for cleaning out numerous files that are not part of the project.
# These remove the following:
# - clean: delete all binary files and temporary output, and the following two.
# - outputclean: delete all *vtk output files
# - dependclean: delete .deps directories and their content

# TODO: rewrite the entire thing ...

# $(context).cleanpattern    := .all.o *.a *.o *.d *.so *.gch *.exe *.exe.stackdump *.out *.out.stackdump OUTPUT_CPPLINT.txt callgrind.out.* 
# $(context).outputcleanpattern := *.vtk
# $(context).depcleanpattern := .deps

# $(context).cleanfiles       := $(patsubst %, $(CURDIR)/%, $(wildcard $($(context).cleanpattern   )) )
# $(context).outputcleanfiles := $(patsubst %, $(CURDIR)/%, $(wildcard $($(context).outputcleanpattern)) )
# $(context).depcleanfiles    := $(patsubst %, $(CURDIR)/%, $(wildcard $($(context).depcleanpattern)) )

# cleanfiles    += $($(context).cleanfiles)
# outputcleanfiles += $($(context).outputcleanfiles)
# depcleanfiles += $($(context).depcleanfiles)

CMD_CLEAN       = rm -f .all.o *.a *.o *.d *.so *.json *.gch OUTPUT_CPPLINT.txt callgrind.out.* *.exe *.exe.stackdump *.out *.out.stackdump 
CMD_OUTPUTCLEAN = rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk ./*.svg ./*/*.svg ./*/*/*.svg ./*.tex ./*/*.tex ./*/*/*.tex
CMD_DEPCLEAN    = if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; rmdir .deps/; fi 

.PHONY: clean vktclean dependclean
.PHONY: $(context).clean $(context).vktclean $(context).dependclean

clean:       $(context).clean
outputclean: $(context).outputclean
dependclean: $(context).dependclean

$(context).clean $(context).outputclean $(context).dependclean: mycontext    := $(context)
$(context).clean $(context).outputclean $(context).dependclean: mycontextdir := $(contextdir)

$(context).clean: 
#	@-echo $(PWD)
	@-cd $(mycontextdir); $(CMD_CLEAN); $(CMD_OUTPUTCLEAN); $(CMD_DEPCLEAN); 

$(context).outputclean: 
#	@-echo $(PWD)
	@-cd $(mycontextdir); $(CMD_OUTPUTCLEAN); 

$(context).dependclean: 
#	@-echo $(PWD)
	@-cd $(mycontextdir); $(CMD_DEPCLEAN); 


