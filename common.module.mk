
# This file contains makefile rules for object creation
# that can be applied in every single source directory
# They create the dependency auxiliary files, 
# the object files, and the shared libraries. 

######################################################################################
# EXPECTED VARIABLES:
# - projectdir : path/to/project/directory
# - moddir     : path/to/module/directory
# - module     : name of the module 
ifndef module
$(error Expect 'module')
endif
ifndef moddir
$(error Expect 'moddir')
endif


######################################################################################
# Set the variables for this file 
# determine whether to use static or dynamic linking 

$(module).depdir  := $(moddir)/.deps

$(module).sources      := $(wildcard $(moddir)/*.cpp)
$(module).headers      := $(wildcard $(moddir)/*.hpp)
$(module).objects      := $(patsubst %.cpp,%.o,$($(module).sources))
$(module).dependencies := $(patsubst %.cpp,$($(module).depdir)/%.d,$(notdir $($(module).sources)))

$(module).libraryobject := $(moddir)/lib$(module).o
$(module).sharedlibrary := $(moddir)/lib$(module).so
$(module).staticlibrary := $(moddir)/lib$(module).a

# Set target specific variables for all recipes
$(module).%: mymodule := $(module)
$(module).%: mymoddir := $(moddir)

###################################################################################################
# All object files that are compiled (and their descendants) also depend on the makefiles 

$(moddir)/*.o $(moddir)/.all.o: $(moddir)/makefile
$(moddir)/*.o $(moddir)/.all.o: $(projectdir)/makefile $(projectdir)/common.compile.mk $(projectdir)/common.module.mk 


###################################################################################################
# How to compile the .o files and .all.o file 

$($(module).depdir):
	@echo $@
	@"mkdir" -p $@

# We generate dependency files during compilation using the following compiler flags 
# -MT $@ : sets the target of the makefile rule, stripped of any directory components
# -MP    : add dependencies as phony targets
# -MF ...: sets the output file for the rules 
# -MMD   : list headers as a by-product of compiling, excluding system headers
# alternatively,
# -MM    : list headers, excluding system headers, and do not compile anything

DEPFLAGS = -MT $@ -MMD -MP -MF $($(mymodule).depdir)/$(notdir $*.d)

$(moddir)/$($(module).objects): mymodule := $(module)
$(moddir)/$($(module).objects): mymoddir := $(moddir)
$(moddir)/$($(module).objects): %.o: %.cpp $($(module).depdir)/%.d | $($(module).depdir)
#	 @echo module        $(mymodule)
#	 @echo mymoddir      $(mymoddir)
#	 @echo module.depdir $(mymodule).depdir
#	 @echo depdir        $($(mymodule).depdir) 
#	 @echo sources       $($(mymodule).sources) 
#	 @echo headers       $($(mymodule).headers) 
#	 @echo objects       $($(mymodule).objects) 
#	 @echo dependencies  $($(mymodule).dependencies) 
#	 @echo sobase        $($(mymodule).sharedlibrarybasename)
#	 @echo libobject     $($(mymodule).libraryobject)
#	 @echo so            $($(mymodule).sharedlibrary)
#	 @echo a             $($(mymodule).staticlibrary)
#	 @echo target:       $@
#	 @echo target(red.)  $*
#	 @echo source file:  $<
#	 @echo all prereq:   $^
#	 @echo DEPFLAGS:     $(DEPFLAGS)
#	 @echo $($(mymodule).depdir)/$(notdir $*.d)
	@echo Compiling and listing dependencies: $@ 
#	$(CXX) $(CXXFLAGS) $(CPPFLAGS)                       $< -c -o $@ $(DEPFLAGS)
	@touch $@

# Set target specific variables for all recipes
# $(module).%: mymodule := $(module)
# $(module).%: mymoddir := $(moddir)

$(moddir)/.all.o: mymoddir := $(module)
$(moddir)/.all.o: mymoddir := $(moddir)

$(moddir)/.all.o: $($(module).sources) $(moddir)/.all.cpp $($(module).depdir)/.all.d | $($(module).depdir)
#	 @echo module        $(mymodule)
#	 @echo mymoddir      $(mymoddir)
#	 @echo module.depdir $(mymodule).depdir
#	 @echo depdir        $($(mymodule).depdir) 
#	 @echo sources       $($(mymodule).sources) 
#	 @echo headers       $($(mymodule).headers) 
#	 @echo objects       $($(mymodule).objects) 
#	 @echo dependencies  $($(mymodule).dependencies) 
#	 @echo sobase        $($(mymodule).sharedlibrarybasename)
#	 @echo libobject     $($(mymodule).libraryobject)
#	 @echo so            $($(mymodule).sharedlibrary)
#	 @echo a             $($(mymodule).staticlibrary)
#	 @echo target:       $@
#	 @echo target_red    $*
#	 @echo source file:  $<
#	 @echo all prereq:   $^
#	 @echo DEPFLAGS:     $(DEPFLAGS)
#	 @echo Compiling and listing dependencies: $($(mymodule).libraryobject)
#	 @echo Compiling and listing dependencies: $@
	@echo Compiling and listing dependencies: $(mymoddir)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(mymoddir)/.all.cpp -c -o $@  $(DEPFLAGS)
	@touch $@

$($(module).depdir)/.all.d:

$($(module).dependencies):

-include $($(module).depdir)/.all.d
-include $($(module).dependencies)


# We can also generate the dependency files without compiling anything 

.PHONY: make_dependencies $(module).make_dependencies
make_dependencies: $(module).make_dependencies

$(module).make_dependencies: $($(module).depdir)
	@for item in $($(mymodule).sources); do \
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $$item -MM -MP -MF $$item.d; \
	done

# How to the .o files into library objects 

$($(module).libraryobject): $(moddir)/.all.o
	@echo Library object: $@
	@cp $(mymoddir)/.all.o $($(mymodule).libraryobject)

$($(module).sharedlibrary): $($(module).libraryobject)
	@echo Shared library: $@
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

$($(module).staticlibrary): $($(module).libraryobject)
	@echo Static library: $@
	@ar crs $($(mymodule).staticlibrary) $($(mymodule).libraryobject)




###################################################################################################
# Finally, determine the build target depending on whether shared libraries are poossible or not

.PHONY: buildobjects buildso builda
.PHONY: $(module).buildobjects $(module).buildso $(module).builda
.PHONY: $(module).build

$(module).buildobjects: $($(module).objects)
$(module).buildso:      $($(module).sharedlibrary)
$(module).builda:       $($(module).staticlibrary)

buildobjects: $(module).buildobjects
buildso:      $(module).buildso
builda:       $(module).builda


ifeq ($(LINKINGTYPE),objectfile)
$(module).build: $($(module).libraryobject)
else
ifeq ($(LINKINGTYPE),static)
$(module).build: $(module).builda
else 
ifeq ($(LINKINGTYPE),dynamic)
$(module).build: $(module).buildso
else 
$(error Unknown linkingtype $(LINKINGTYPE))
endif 
endif 
endif

# ifeq ($(OS),Windows_NT)
# $(module).build: $(module).builda
# else
# $(module).build: $(module).builda $(module).buildso
# endif





##########################################################################################
# List several objects. Read-only

.PHONY:  list_of_objects $(module).list_of_objects
.SILENT: list_of_objects $(module).list_of_objects

list_of_objects: $(module).list_of_objects

$(module).list_of_objects: 
	@echo Project directory: $(projecdir);
	@echo Module name:       $(mymodule);
	@echo Module directory:  $(mymoddir);
	@echo Sources:           $($(mymodule).sources);
	@echo Headers:           $($(mymodule).headers);
	@echo Static lib:        $($(mymodule).staticlibrary);
	@echo Shared lib:        $($(mymodule).sharedlibrary);
	@echo Objects:           $($(mymodule).objects);
	@echo Dependencies:      $($(mymodule).dependencies);
	@echo Build:             $($(mymodule).build);


ifneq ($(module),external)
##########################################################################################
# Check whether the header files have correct syntax. Read-only.

$(module).headerchecks := $(patsubst %.hpp,check-%.hpp,$($(module).headers))

.PHONY: checkheaders $($(module).headerchecks)
checkheaders: $($(module).headerchecks)
$($(module).headerchecks): check-%.hpp : 
	$(info Check header: $*.hpp)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $*.hpp -fsyntax-only




##########################################################################################
# Check whether the source files have correct syntax. Read-only.

$(module).sourcechecks := $(patsubst %.cpp,check-%.cpp,$($(module).sources))

.PHONY: checksources $($(module).sourcechecks)
checksources: $($(module).sourcechecks)
$($(module).sourcechecks): check-%.cpp : 
	$(info Check source: $*.cpp)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $*.cpp -fsyntax-only



##########################################################################################
# Target 'check' is a generic test. Currently, it defaults to 'tidy'

.PHONY: check
check: tidy

##########################################################################################
# Apply clang-tidy to all cpp and hpp files in the directory. Read-only.

.PHONY: tidy $(module).tidy
tidy: $(module).tidy
$(module).tidy:
	clang-tidy $($(mymodule).sources) --config-file=$(projectdir)/.Tools/clang-tidy.yaml -- -std=c++17 # -fno-exceptions


##########################################################################################
# Apply cppcheck to all cpp and hpp files in the directory. Read-only. 

.PHONY: cppcheck $(module).cppcheck
cppcheck: $(module).cppcheck
$(module).cppcheck:
	cppcheck -i ./.playground/ -i ./.legacy/ -i ./external/ \
	--enable=warning,style,performance,portability \
	--suppress=duplicateCondition --suppress=toomanyconfigs --suppress=sizeofFunctionCall --suppress=variableScope --suppress=missingReturn --suppress=assertWithSideEffect --suppress=useStlAlgorithm --suppress=knownConditionTrueFalse --suppress=unsignedPositive \
	--std=c++17 -q $(mymoddir)/*pp


##########################################################################################
# Regex several useful things. Read-only. 
# - find trailing white spaces 
# - find non-ASCII characters 
# - find consecutive spaces 

.PHONY: grepissues $(module).grepissues
grepissues: $(module).grepissues
# $(module).grepissues: mymodule := $(module)
# $(module).grepissues: mymoddir := $(moddir)
$(module).grepissues:
	@echo Search trailing whitespace...
	@-grep --line-number --color '\s+$$' -r $(mymoddir)/*pp
#	@echo Search non-ASCII characters...
#	@-grep --line-number --color '[^\x00-\x7F]' -r $(mymoddir)/*pp
#	@echo Find consecutive spaces...
#	@-grep --line-number --color '\b\s{2,}' -r $(mymoddir)/*pp
#	@echo Find standard asserts...
#	@-grep --line-number --color 'assert(' $(mymoddir)/*pp
#	@echo Find usage of 'cout' ...
#	@-grep --line-number --color 'cout' $(mymoddir)/*pp
#	@echo Find floating-point numbers ...
#	@-grep --line-number --color -E '\.*[0-9]' $(mymoddir)/*pp
#	@-grep --line-number --color -E '(0-9)e' $(mymoddir)/*pp
#	@-grep --line-number --color -E '([0-9]+e[0-9]+)|([0-9]+\.[0-9]+)|((+-\ )\.[0-9]+)|((+-\ )[0-9]+\.)' $(mymoddir)/*pp

endif



##########################################################################################
# Commands for cleaning out numerous files that are not part of the project.
# These remove the following:
# - clean:       delete all binary files and temporary output, and the following two.
# - outputclean: delete all *vtk output files
# - dependclean: delete .deps directories and their content

CMD_CLEAN       = rm -f .all.o .all.json .json *.a *.o *.d *.so *.json *.gch OUTPUT_CPPLINT.txt callgrind.out.* *.exe *.exe.stackdump *.out *.out.stackdump 
CMD_OUTPUTCLEAN = rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk ./*.svg ./*/*.svg ./*/*/*.svg
CMD_DEPCLEAN    = if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; rmdir .deps/; fi 

.PHONY: clean outputclean dependclean
.PHONY: $(module).clean $(module).outputclean $(module).dependclean

clean:       $(module).clean
outputclean: $(module).outputclean
dependclean: $(module).dependclean

$(module).clean $(module).outputclean $(module).dependclean: mymodule := $(module)
$(module).clean $(module).outputclean $(module).dependclean: mymoddir := $(moddir)

$(module).clean: 
#	@-echo $(PWD)
	@-cd $(mymoddir); $(CMD_CLEAN); $(CMD_OUTPUTCLEAN); $(CMD_DEPCLEAN); 

$(module).outputclean: 
#	@-echo $(PWD)
	@-cd $(mymoddir); $(CMD_OUTPUTCLEAN); 

$(module).dependclean: 
#	@-echo $(PWD)
	@-cd $(mymoddir); $(CMD_DEPCLEAN); 


