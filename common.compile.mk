
# 
# This file contains definitions of variables 
# to be used in the compilation process 
# for objects and shared libraries 
# 
# This file produces the variables 
# 	- CXX
# 	- CXXFLAGS
# 	- CPPFLAGS
# 	- LDFLAGS
# 	- LDLIBS
# 
# Note that these include several preprocessor variables
# Furthermore, it sets:
# 
#   - LINKINGTYPE
# 
#  

################################################################################
# EXPECTED VARIABLES:
# - projectdir : path/to/project/directory
ifndef projectdir
$(error Expect 'projectdir')
endif





# Uncomment your choice of compiler below
FLAG_CXX := CLANG
# FLAG_CXX := GCC
# FLAG_CXX := ICC


# Do you want to DISABLE the general assert macro?
# Uncomment the following line to disable the general assert macro
# FLAG_DISABLE_ASSERTIONS=yes

# Do you want the standard library assert macro instead of the custom one?
# Uncomment the following line to use the standard library assert macro 
# FLAG_USE_ORIGINAL_ASSERT_MACRO=yes

# Do you want assert messages to be discarded?
# Uncomment the following line to simplify the debugging macros 
# FLAG_DISCARD_ASSERT_MESSAGES=yes

# Do you want to DISABLE checking of meshes?
# Uncomment the following line to disable extensive check routines for meshes
FLAG_DISABLE_CHECK_MESHES=yes

# Do you want to ENABLE the standard library debugging flags 
# Uncomment the following line to enable the standard library debugging flags 
# FLAG_DISABLE_STDLIBDEBUG=yes

# Do you want to DISABLE the custom logging framework
# in favor of standard library routines?
# Uncomment the following line for that
# FLAG_USE_PRIMITIVE_LOGGING=yes

# Do you want to DISABLE excpetion handling?
# Uncomment the following line to disable exception handling
FLAG_NO_EXCEPTIONS=yes

# Do you want to compile with all optimization flags enabled?
# Uncomment the following line to have this done so
# FLAG_DO_OPTIMIZE=yes

# Do you want to ENABLE the use of openMP?
# Uncomment the following line to enable compilation with openMP
FLAG_ENABLE_OPENMP=yes

# Do you want to ENABLE excessive warning options?
# Uncomment the following line to enable excessive warning options
FLAG_EXCESSIVE_WARNINGS=yes

# Do you want to ENABLE either extended precision or single precision?
# Uncomment one the following lines to switch from double precision
# to either extended precision or single precision
# FLAG_DO_USE_EXTENDED_PRECISION=yes
# FLAG_DO_USE_SINGLE_PRECISION=yes

# Logging output in color?
FLAG_COLORED_OUTPUT=yes

# Do you want to ENABLE the Clang sanitizer?
# Uncomment the following line to enable compilation with the Clang sanitizer
# FLAG_DO_USE_SANITIZER=yes

# Do you want to enable static analysis during the compilation process
# Uncomment the following line to enable static analysis
# FLAG_DO_STATICANALYSIS=yes

# Do you want to ENABLE the use of tcmalloc?
# Uncomment the following line to enable tcmalloc
# FLAG_USE_TCMALLOC=yes

# Do you want to DISABLE embedding of Debug information?
# Uncomment the following line to have no debug information included
# FLAG_NO_DEBUGINFO=yes

# Do you want to ENABLE the backtracer for debugging?
# Uncomment the following line to enable backtracing
# FLAG_USE_BACKTRACER=yes

# Do you want to strip unused symbols from the executables?
# Uncomment the following line to accomplish this
# FLAG_DO_STRIP=yes

# Do you want to ENABLE profile generation? 
# Uncomment the following line to enable profile generation at every run.
# FLAG_DO_PROFILE=yes

# Choose the linker by uncommenting, or leave commented for default linker 
# LDFLAGS += -fuse-ld=bfd
# LDFLAGS += -fuse-ld=lld
# LDFLAGS += -fuse-ld=gold
# LDFLAGS += -fuse-ld=mold





###############################################
#                                             #
#       Post-process the option choices       #
#                                             #
###############################################

# Use this file to overwrite the default settings above on a local machine
# At this point, only the flags above will be set.
-include $(projectdir)/OVERWRITE.COMPILE.mk

# If we are in RELEASE_MODE then set the following flags 

ifdef RELEASE_MODE
FLAG_DISABLE_CHECK_MESHES=yes
FLAG_DISABLE_STDLIBDEBUG=yes
FLAG_DISABLE_ASSERTIONS=yes
FLAG_DO_OPTIMIZE=yes
FLAG_ENABLE_OPENMP=yes
FLAG_NO_DEBUGINFO=yes
FLAG_NO_EXCEPTIONS=yes
FLAG_DO_STRIP=yes
FLAG_USE_PRIMITIVE_LOGGING=yes
endif

ifdef FLAG_DO_USE_EXTENDED_PRECISION
ifdef FLAG_DO_USE_SINGLE_PRECISION
$(error Extended and single precision requested at the same time)
endif
endif

ifneq ($(FLAG_CXX),GCC)
ifdef FLAG_DO_USE_SINGLE_PRECISION
$(error Single precision floating-point setup only with GCC)
endif
endif





###############################################
#                                             #
#         Set the compiler command            #
#       (see also language std below)         #
###############################################

ifeq ($(FLAG_CXX),GCC)

  CXX := g++ -D__USE_MINGW_ANSI_STDIO=1
  #-ftime-report
  #-fuse-ld=lld
  
else ifeq ($(FLAG_CXX),CLANG)

  CXX := clang++ 
  # -ftime-trace
  # -stdlib=libstdc++ 

else ifeq ($(FLAG_CXX),ICC)

  CXX := icc 

else

  $(error No compiler recognized)

endif

###############################################
#                                             #
#           Language std settings             #
#                                             #
###############################################

CXXFLAGS_LANG := -std=c++20 -pedantic -fno-rtti -D_LIBCPP_REMOVE_TRANSITIVE_INCLUDES 












###############################################
#                                             #
#               Optimization                  #
#                                             #
###############################################

# If optimization is enabled, then set a number of flags 
# In the absence of optimization, we set O0

CXXFLAGS_OPTIMIZE:=

# CXXFLAGS_OPTIMIZE += -march=native -mtune=native 

ifeq ($(FLAG_DO_OPTIMIZE),yes)

	ifeq ($(FLAG_CXX),ICC)
		CXXFLAGS_OPTIMIZE += -march=core-avx2
		CXXFLAGS_OPTIMIZE += -intel-optimized-headers 
		CXXFLAGS_OPTIMIZE += -xHOST -O3 -Ofast -ipo -no-prec-div -fp-model fast=2
		CXXFLAGS_OPTIMIZE += 
		CXXFLAGS_OPTIMIZE += 
		CXXFLAGS_OPTIMIZE += 
	else 
# wierd warnings appear at LTO and O1+ ....
		CXXFLAGS_OPTIMIZE += -flto
		CXXFLAGS_OPTIMIZE += -Ofast  
		CXXFLAGS_OPTIMIZE += -fshort-enums
		ifeq ($(FLAG_CXX),GCC)
			CXXFLAGS_OPTIMIZE += -fno-fat-lto-objects
# 			CXXFLAGS_OPTIMIZE += -finline-limit=1200
# 			CXXFLAGS_OPTIMIZE += -fno-signed-zeros -fno-trapping-math -fassociative-math
			CXXFLAGS_OPTIMIZE += -fmerge-all-constants
			CXXFLAGS_OPTIMIZE += -fipa-pta 
			CXXFLAGS_OPTIMIZE += -fdevirtualize-speculatively
# Loop unrolling is very speculative			
			CXXFLAGS_OPTIMIZE += -funroll-loops -fvariable-expansion-in-unroller -floop-nest-optimize
			CXXFLAGS_OPTIMIZE += -fomit-frame-pointer
#			CXXFLAGS_OPTIMIZE += -malign-double 
		endif
		ifeq ($(FLAG_CXX),CLANG)
#			 CXXFLAGS_OPTIMIZE += -finline-limit=1200
		endif
	endif
else
	CXXFLAGS_OPTIMIZE += -O0
endif


# Do we apply OpenMP?

ifeq ($(FLAG_ENABLE_OPENMP),yes)
	CXXFLAGS_OPTIMIZE += -fopenmp
endif


# Do we strip debug information?

ifeq ($(FLAG_CXX),GCC) 
ifeq ($(FLAG_DO_STRIP),yes)
	CXXFLAGS_OPTIMIZE += -ffunction-sections -fdata-sections -Wl,--gc-sections -Wl,--strip-all 
endif
endif




###############################################
#                                             #
#        Misc Code generation options         #
#                                             #
###############################################

CXXFLAGS_CODEGEN := 

ifeq ($(FLAG_NO_EXCEPTIONS),yes)
	CXXFLAGS_CODEGEN += -fno-exceptions
endif

ifeq ($(FLAG_DO_USE_SINGLE_PRECISION),yes)
	CXXFLAGS_CODEGEN += -fsingle-precision-constant
endif

# If we are NOT on Windows, then use position-independent code. Also, avoid the procedure linkage table 
ifneq ($(OS),Windows_NT)
	CXXFLAGS_CODEGEN += -fpic -fno-plt 
endif

# CXXFLAGS_CODEGEN += -fvisibility=default

CXXFLAGS_CODEGEN += -fvisibility-inlines-hidden

# enums can only take values in the enumeration type 
CXXFLAGS_CODEGEN += -fstrict-enums

# CXXFLAGS_CODEGEN += -fstrict-eval-order





















###############################################
#                                             #
#        Compiler warning settings            #
#                                             #
###############################################

CXXFLAGS_WARNINGS := 
CXXFLAGS_WARNINGS += -Wall -Wextra -Wpedantic 

ifeq ($(FLAG_EXCESSIVE_WARNINGS),yes)

	CXXFLAGS_WARNINGS += -Wodr  
	
	CXXFLAGS_WARNINGS += -Wredundant-decls
	CXXFLAGS_WARNINGS += -Wmissing-declarations
	CXXFLAGS_WARNINGS += -Wmissing-field-initializers
	CXXFLAGS_WARNINGS += -Wctor-dtor-privacy
	CXXFLAGS_WARNINGS += -Wnon-virtual-dtor
	CXXFLAGS_WARNINGS += -Woverloaded-virtual
	
	CXXFLAGS_WARNINGS += -Wmisleading-indentation
	
	CXXFLAGS_WARNINGS += -Wsign-promo
	CXXFLAGS_WARNINGS += -Wcast-align
	CXXFLAGS_WARNINGS += -Wcast-qual
	# # DISABLED: CXXFLAGS_WARNINGS += -Wold-style-cast
	
	CXXFLAGS_WARNINGS += -Wundef
	CXXFLAGS_WARNINGS += -Wunused 
	
	CXXFLAGS_WARNINGS += -Wformat=2
	CXXFLAGS_WARNINGS += -Wformat-nonliteral
	CXXFLAGS_WARNINGS += -Wformat-security
	CXXFLAGS_WARNINGS += -Wformat-y2k
	
	# # DISABLED: CXXFLAGS_WARNINGS += -Wcomma-subscript
	# # DISABLED: CXXFLAGS_WARNINGS += -Wlifetime
	# # DISABLED: CXXFLAGS_WARNINGS += -Weffc++

	ifeq ($(FLAG_CXX),GCC)

		# # DISABLED: CXXFLAGS_WARNINGS += -Wabi
		CXXFLAGS_WARNINGS += -Waggressive-loop-optimizations
		CXXFLAGS_WARNINGS += -Walloca
		# # DISABLED: CXXFLAGS_WARNINGS += -Walloc-size
		CXXFLAGS_WARNINGS += -Walloc-zero
		CXXFLAGS_WARNINGS += -Warith-conversion
		CXXFLAGS_WARNINGS += -Warray-bounds
		CXXFLAGS_WARNINGS += -Warray-bounds=2
		CXXFLAGS_WARNINGS += -Wattribute-alias=2
		# # unknown: CXXFLAGS_WARNINGS += -Wbidi-chars=any
		# # unknown: CXXFLAGS_WARNINGS += -Wcalloc-transposed-args
		CXXFLAGS_WARNINGS += -Wcast-qual
		CXXFLAGS_WARNINGS += -Wcast-align
		CXXFLAGS_WARNINGS += -Wcast-align=strict 
		CXXFLAGS_WARNINGS += -Wconversion
		CXXFLAGS_WARNINGS += -Wdangling-else 
		CXXFLAGS_WARNINGS += -Wdate-time
		CXXFLAGS_WARNINGS += -Wdisabled-optimization
		CXXFLAGS_WARNINGS += -Wduplicated-branches 
		CXXFLAGS_WARNINGS += -Wduplicated-cond
		CXXFLAGS_WARNINGS += -Wexpansion-to-defined
		CXXFLAGS_WARNINGS += -Wformat=2
		CXXFLAGS_WARNINGS += -Wformat-overflow=2
		CXXFLAGS_WARNINGS += -Wformat-signedness
		CXXFLAGS_WARNINGS += -Wformat-truncation=1
		# # disabled: CXXFLAGS_WARNINGS += -Wformat-truncation=2
		CXXFLAGS_WARNINGS += -Wfloat-conversion
		CXXFLAGS_WARNINGS += -Wfloat-equal
		CXXFLAGS_WARNINGS += -Wfree-nonheap-object
		CXXFLAGS_WARNINGS += -Wimplicit-fallthrough=4
		# # DISABLED: CXXFLAGS_WARNINGS += -Winfinite-recursion
		CXXFLAGS_WARNINGS += -Winit-self
		# # DISABLED: CXXFLAGS_WARNINGS += -Winline 
		CXXFLAGS_WARNINGS += -Winline 
		CXXFLAGS_WARNINGS += -Wint-in-bool-context
		CXXFLAGS_WARNINGS += -Winvalid-pch
		# # unknown: CXXFLAGS_WARNINGS += -Winvalid-utf8
		CXXFLAGS_WARNINGS += -Wlogical-op 
		# # DISABLED: CXXFLAGS_WARNINGS += -Wlong-long
		CXXFLAGS_WARNINGS += -Wmisleading-indentation 
		CXXFLAGS_WARNINGS += -Wmissing-declarations
		CXXFLAGS_WARNINGS += -Wmissing-include-dirs
		CXXFLAGS_WARNINGS += -Wmultiple-inheritance 
		# # DISABLED: CXXFLAGS_WARNINGS += -Wnrvo
		CXXFLAGS_WARNINGS += -Wnull-dereference
		# # DISABLED: CXXFLAGS_WARNINGS += -Wpadded
		CXXFLAGS_WARNINGS += -Wparentheses
		CXXFLAGS_WARNINGS += -Wpointer-arith
		CXXFLAGS_WARNINGS += -Wpointer-arith
		CXXFLAGS_WARNINGS += -Wredundant-decls
		CXXFLAGS_WARNINGS += -Wreturn-local-addr
		# # DISABLED: CXXFLAGS_WARNINGS += -Wreturn-mismatch
		# # DISABLED: CXXFLAGS_WARNINGS += -Wshadow
		# # DISABLED: CXXFLAGS_WARNINGS += -Wshadow=global
		CXXFLAGS_WARNINGS += -Wshift-negative-value
		CXXFLAGS_WARNINGS += -Wshift-overflow=2
		CXXFLAGS_WARNINGS += -Wsign-compare
		CXXFLAGS_WARNINGS += -Wsign-conversion
		CXXFLAGS_WARNINGS += -Wsign-promo 
		CXXFLAGS_WARNINGS += -Wstrict-aliasing=3
		# # DISABLED: CXXFLAGS_WARNINGS += -Wstrict-flex-arrays
		# # DISABLED: CXXFLAGS_WARNINGS += -Wstrict-overflow=5
		CXXFLAGS_WARNINGS += -Wstringop-overflow=4
		# # DISABLED: CXXFLAGS_WARNINGS += -Wsuggest-attribute
		CXXFLAGS_WARNINGS += -Wsuggest-final-types 
		CXXFLAGS_WARNINGS += -Wsuggest-final-methods 
		CXXFLAGS_WARNINGS += -Wsuggest-override
		CXXFLAGS_WARNINGS += -Wswitch-default 
		CXXFLAGS_WARNINGS += -Wswitch-enum
		CXXFLAGS_WARNINGS += -Wsync-nand
		CXXFLAGS_WARNINGS += -Wtautological-compare
		# # DISABELD: CXXFLAGS_WARNINGS += -Wtraditional
		CXXFLAGS_WARNINGS += -Wtrampolines
		# # DISABELD: CXXFLAGS_WARNINGS += -Wtrivial-auto-var-init
		CXXFLAGS_WARNINGS += -Wuninitialized
		CXXFLAGS_WARNINGS += -Wunused 
		CXXFLAGS_WARNINGS += -Wunused-const-variable=1
		# # DISABLED: CXXFLAGS_WARNINGS += -Wunused-macros
		# # DISABLED: CXXFLAGS_WARNINGS += -Wuseless-cast
		CXXFLAGS_WARNINGS += -Wunsafe-loop-optimizations
		CXXFLAGS_WARNINGS += -Wundef
		CXXFLAGS_WARNINGS += -Wvector-operation-performance
		CXXFLAGS_WARNINGS += -Wvirtual-inheritance
		CXXFLAGS_WARNINGS += -Wwrite-strings
		CXXFLAGS_WARNINGS += -Wzero-length-bounds
		CXXFLAGS_WARNINGS += -Wzero-as-null-pointer-constant
		
		# # check which options there are ... 
		# # TODO What is stack smashing???? HSA????
		# # TODO Read the format warnings 

		CXXFLAGS_WARNINGS += 
		
	else ifeq ($(FLAG_CXX),CLANG)

# 		CXXFLAGS_WARNINGS += -Wdangling-reference
		
		CXXFLAGS_WARNINGS += -Wabstract-vbase-init
		CXXFLAGS_WARNINGS += -Walloca
		CXXFLAGS_WARNINGS += -Wanon-enum-enum-conversion
		CXXFLAGS_WARNINGS += -Werror-implicit-function-declaration 
		CXXFLAGS_WARNINGS += -Wabsolute-value 
		CXXFLAGS_WARNINGS += -Wanon-enum-enum-conversion 
		CXXFLAGS_WARNINGS += -Warray-bounds-pointer-arithmetic
		CXXFLAGS_WARNINGS += -Wasm-operand-widths
		CXXFLAGS_WARNINGS += -Wassign-enum
		CXXFLAGS_WARNINGS += -Watomic-implicit-seq-cst
		CXXFLAGS_WARNINGS += -Watomic-properties
		# # unknown: CXXFLAGS_WARNINGS += -Wauto-decl-extensions
		CXXFLAGS_WARNINGS += -Wbad-function-cast
		# # disabled: CXXFLAGS_WARNINGS += -Wbinary-literal
		CXXFLAGS_WARNINGS += -Wbind-to-temporary-copy
		CXXFLAGS_WARNINGS += -Wbit-int-extension
		CXXFLAGS_WARNINGS += -Wbitfield-enum-conversion
		CXXFLAGS_WARNINGS += -Wbitwise-instead-of-logical
		CXXFLAGS_WARNINGS += -Wbitwise-op-parentheses
		CXXFLAGS_WARNINGS += -Wbool-operation
		# # disabled: CXXFLAGS_WARNINGS += -Wc++-compat
		CXXFLAGS_WARNINGS += -Wc++11-narrowing
		# -Wc++11-compat-pedantic
		# -Wc++11-compat-reserved-user-defined-literal
		# -Wc++11-extensions
		# -Wc++11-extra-semi
		# -Wc++11-long-long
		# -Wc++11-narrowing
		# -Wc++11-narrowing-const-reference
		# # Skipping language version specific warnings ...
		CXXFLAGS_WARNINGS += -Wcall-to-pure-virtual-from-ctor-dtor
		CXXFLAGS_WARNINGS += -Wcalled-once-parameter
		CXXFLAGS_WARNINGS += -Wcast-align
		CXXFLAGS_WARNINGS += -Wcast-function-type
		# # unknown: CXXFLAGS_WARNINGS += -Wcast-function-type-strict
		CXXFLAGS_WARNINGS += -Wcast-qual
		CXXFLAGS_WARNINGS += -Wchar-subscripts
		CXXFLAGS_WARNINGS += -Wclass-varargs
		CXXFLAGS_WARNINGS += -Wcomma
		CXXFLAGS_WARNINGS += -Wcomment
		CXXFLAGS_WARNINGS += -Wcomplex-component-init
		CXXFLAGS_WARNINGS += -Wcompound-token-split
		CXXFLAGS_WARNINGS += -Wconditional-type-mismatch
		CXXFLAGS_WARNINGS += -Wconditional-uninitialized
		CXXFLAGS_WARNINGS += -Wconsumed
		CXXFLAGS_WARNINGS += -Wconversion 
		# # disable: CXXFLAGS_WARNINGS += -Wcovered-switch-default
		CXXFLAGS_WARNINGS += -Wctad-maybe-unsupported
		CXXFLAGS_WARNINGS += -Wcuda-compat
		CXXFLAGS_WARNINGS += -Wdate-time
		CXXFLAGS_WARNINGS += -Wdelegating-ctor-cycles
		CXXFLAGS_WARNINGS += -Wdelete-non-abstract-non-virtual-dtor
		CXXFLAGS_WARNINGS += -Wdelete-non-virtual-dtor
		CXXFLAGS_WARNINGS += -Wdelimited-escape-sequence-extension
		CXXFLAGS_WARNINGS += -Wdeprecated
		CXXFLAGS_WARNINGS += -Wdirect-ivar-access
		CXXFLAGS_WARNINGS += -Wdisabled-macro-expansion
		CXXFLAGS_WARNINGS += -Wdollar-in-identifier-extension
		# # DISABLED: CXXFLAGS_WARNINGS += -Wdouble-promotion #disabled
		#	CXXFLAGS_WARNINGS += -Wdtor-name  TODO: Is this one actually defined???
		CXXFLAGS_WARNINGS += -Wduplicate-decl-specifier 
		CXXFLAGS_WARNINGS += -Wduplicate-enum 
		CXXFLAGS_WARNINGS += -Wduplicate-method-arg 
		CXXFLAGS_WARNINGS += -Wduplicate-method-match 
		CXXFLAGS_WARNINGS += -Wembedded-directive
		CXXFLAGS_WARNINGS += -Wempty-init-stmt
		CXXFLAGS_WARNINGS += -Wempty-translation-unit
		CXXFLAGS_WARNINGS += -Wenum-compare-conditional
		CXXFLAGS_WARNINGS += -Wexceptions
		# # disabled: CXXFLAGS_WARNINGS += -Wexit-time-destructors
		CXXFLAGS_WARNINGS += -Wexpansion-to-defined
		CXXFLAGS_WARNINGS += -Wextra-semi
		# # disabled: CXXFLAGS_WARNINGS += -Wextra-semi-stmt
		
		CXXFLAGS_WARNINGS += -Wfixed-enum-extension
		CXXFLAGS_WARNINGS += -Wflexible-array-extensions
		# # DISABLED: CXXFLAGS_WARNINGS += -Wfloat-equal
		CXXFLAGS_WARNINGS += -Wfloat-overflow-conversion
		CXXFLAGS_WARNINGS += -Wfloat-conversion
		CXXFLAGS_WARNINGS += -Wfloat-zero-conversion
		CXXFLAGS_WARNINGS += -Wfor-loop-analysis
		CXXFLAGS_WARNINGS += -Wformat
		CXXFLAGS_WARNINGS += -Wformat=2
		CXXFLAGS_WARNINGS += -Wfortify-source
		CXXFLAGS_WARNINGS += -Wfour-char-constants
		CXXFLAGS_WARNINGS += -Wframe-address
		# # DISABLED: -Wframe-larger-than
		CXXFLAGS_WARNINGS += -Wfuse-ld-path
		CXXFLAGS_WARNINGS += -Wgcc-compat
		# # unknown: CXXFLAGS_WARNINGS += -Wgeneric-type-extension
		# # disabled: CXXFLAGS_WARNINGS += -Wglobal-constructors
		CXXFLAGS_WARNINGS += -Wgnu
		CXXFLAGS_WARNINGS += -Wheader-guard
		CXXFLAGS_WARNINGS += -Wheader-hygiene
		CXXFLAGS_WARNINGS += -Widiomatic-parentheses
		CXXFLAGS_WARNINGS += -Wimplicit
		CXXFLAGS_WARNINGS += -Wimplicit-atomic-properties
		CXXFLAGS_WARNINGS += -Wimplicit-fallthrough
		CXXFLAGS_WARNINGS += -Wimplicit-fallthrough-per-function
		CXXFLAGS_WARNINGS += -Wimplicit-int-conversion
		CXXFLAGS_WARNINGS += -Wimplicit-int-float-conversion
		CXXFLAGS_WARNINGS += -Wimplicit-retain-self
		CXXFLAGS_WARNINGS += -Wimplicit-float-conversion
		CXXFLAGS_WARNINGS += -Wimplicit-function-declaration
		CXXFLAGS_WARNINGS += -Wimport-preprocessor-directive-pedantic
		# # unknown: CXXFLAGS_WARNINGS += -Wincompatible-function-pointer-types-strict
		CXXFLAGS_WARNINGS += -Wincomplete-module
		# # disabled: CXXFLAGS_WARNINGS += -Winconsistent-missing-destructor-override
		CXXFLAGS_WARNINGS += -Winfinite-recursion
		CXXFLAGS_WARNINGS += -Winit-self
		CXXFLAGS_WARNINGS += -Winitializer-overrides
		CXXFLAGS_WARNINGS += -Wint-conversion
		CXXFLAGS_WARNINGS += -Wint-in-bool-context
		CXXFLAGS_WARNINGS += -Winvalid-or-nonexistent-directory
		CXXFLAGS_WARNINGS += -Winvalid-partial-specialization
		# # unknown: CXXFLAGS_WARNINGS += -Winvalid-utf8
		CXXFLAGS_WARNINGS += -Wkeyword-macro
		CXXFLAGS_WARNINGS += -Wlanguage-extension-token
		# # disabled: CXXFLAGS_WARNINGS += -Wlocal-type-template-args
		CXXFLAGS_WARNINGS += -Wlogical-op-parentheses
		CXXFLAGS_WARNINGS += -Wloop-analysis
		#CXXFLAGS_WARNINGS += -Wmain
		CXXFLAGS_WARNINGS += -Wmethod-signatures
		CXXFLAGS_WARNINGS += -Wmicrosoft
		CXXFLAGS_WARNINGS += -Wmisleading-indentation
		CXXFLAGS_WARNINGS += -Wmismatched-tags
		CXXFLAGS_WARNINGS += -Wmissing-braces
		CXXFLAGS_WARNINGS += -Wmissing-field-initializers
		CXXFLAGS_WARNINGS += -Wmissing-method-return-type
		CXXFLAGS_WARNINGS += -Wmissing-noreturn
		CXXFLAGS_WARNINGS += -Wmissing-prototypes
		CXXFLAGS_WARNINGS += -Wmissing-variable-declarations
		# no module warnings and remarks yet ...
		CXXFLAGS_WARNINGS += -Wmost 
		CXXFLAGS_WARNINGS += -Wmove 
		CXXFLAGS_WARNINGS += -Wnarrowing
		CXXFLAGS_WARNINGS += -Wnested-anon-types
		CXXFLAGS_WARNINGS += -Wnewline-eof
		CXXFLAGS_WARNINGS += -Wnon-gcc
		CXXFLAGS_WARNINGS += -Wnon-virtual-dtor
		CXXFLAGS_WARNINGS += -Wnonportable-system-include-path
		CXXFLAGS_WARNINGS += -Wnull-pointer-arithmetic 
		CXXFLAGS_WARNINGS += -Wnull-pointer-subtraction
		CXXFLAGS_WARNINGS += -Wnullability-extension
		CXXFLAGS_WARNINGS += -Wnullable-to-nonnull-conversion
		# Skipping object C warnings ...
		# # disabled: CXXFLAGS_WARNINGS += -Wold-style-cast
		CXXFLAGS_WARNINGS += -Wopenmp
		CXXFLAGS_WARNINGS += -Woption-ignored
		CXXFLAGS_WARNINGS += -Wover-aligned
		CXXFLAGS_WARNINGS += -Woverlength-strings
		CXXFLAGS_WARNINGS += -Woverloaded-virtual
		CXXFLAGS_WARNINGS += -Woverriding-method-mismatch
		
		
		CXXFLAGS_WARNINGS += -Wpacked
		# # unknown: CXXFLAGS_WARNINGS += -Wpacked-non-pod
		# # disabled: CXXFLAGS_WARNINGS += -Wpadded
		CXXFLAGS_WARNINGS += -Wparentheses
		CXXFLAGS_WARNINGS += -Wpessimizing-move
		CXXFLAGS_WARNINGS += -Wpointer-arith
		CXXFLAGS_WARNINGS += -Wpoison-system-directories
		CXXFLAGS_WARNINGS += -Wpragma-pack-suspicious-include
		# Skipping pre comp warnings 
		CXXFLAGS_WARNINGS += -Wprofile-instr-missing
		CXXFLAGS_WARNINGS += -Wquoted-include-in-framework-header
		CXXFLAGS_WARNINGS += -Wrange-loop-analysis
		CXXFLAGS_WARNINGS += -Wredundant-move
		
		CXXFLAGS_WARNINGS += -Wredundant-parens
		CXXFLAGS_WARNINGS += -Wreorder-ctor
		CXXFLAGS_WARNINGS += -Wreserved-identifier
		CXXFLAGS_WARNINGS += -Wreserved-id-macro
		CXXFLAGS_WARNINGS += -Wreserved-user-defined-literal
		CXXFLAGS_WARNINGS += -Wretained-language-linkage
		CXXFLAGS_WARNINGS += -Wreturn-std-move
		
		CXXFLAGS_WARNINGS += -Wself-assign 
		CXXFLAGS_WARNINGS += -Wself-move
		CXXFLAGS_WARNINGS += -Wsemicolon-before-method-body
		# # disabled: CXXFLAGS_WARNINGS += -Wshadow
		# # disabled: CXXFLAGS_WARNINGS += -Wshadow-all
		CXXFLAGS_WARNINGS += -Wshift-sign-overflow
		CXXFLAGS_WARNINGS += -Wshorten-64-to-32
		CXXFLAGS_WARNINGS += -Wsign-compare
		CXXFLAGS_WARNINGS += -Wsign-conversion
		CXXFLAGS_WARNINGS += -Wsigned-enum-bitfield
		# # unknown: CXXFLAGS_WARNINGS += -Wsingle-bit-bitfield-constant-conversion
		CXXFLAGS_WARNINGS += -Wsizeof-array-argument
		CXXFLAGS_WARNINGS += -Wsizeof-array-decay
		CXXFLAGS_WARNINGS += -Wsizeof-array-div
		CXXFLAGS_WARNINGS += -Wsizeof-pointer-div
		CXXFLAGS_WARNINGS += -Wsizeof-pointer-memaccess


		CXXFLAGS_WARNINGS += -Wsometimes-uninitialized
		# # unknown: CXXFLAGS_WARNINGS += -Wsource-uses-openacc
		CXXFLAGS_WARNINGS += -Wsource-uses-openmp
		CXXFLAGS_WARNINGS += -Wstack-exhausted
		CXXFLAGS_WARNINGS += -Wstatic-in-inline
		CXXFLAGS_WARNINGS += -Wstrict-prototypes
		CXXFLAGS_WARNINGS += -Wstring-compare
		CXXFLAGS_WARNINGS += -Wstring-concatenation
		CXXFLAGS_WARNINGS += -Wstring-conversion
		# # unknown: CXXFLAGS_WARNINGS += -Wsuggest-destructor-override
		CXXFLAGS_WARNINGS += -Wsuggest-override
		CXXFLAGS_WARNINGS += -Wsuper-class-method-mismatch
		CXXFLAGS_WARNINGS += -Wswitch-default
		CXXFLAGS_WARNINGS += -Wswitch-enum
		CXXFLAGS_WARNINGS += -Wtautological-bitwise-compare
		CXXFLAGS_WARNINGS += -Wtautological-compare
		CXXFLAGS_WARNINGS += -Wtautological-constant-compare
		CXXFLAGS_WARNINGS += -Wtautological-constant-in-range-compare
		CXXFLAGS_WARNINGS += -Wtautological-constant-out-of-range-compare
		# # unknown: CXXFLAGS_WARNINGS += -Wtautological-negation-compare
		CXXFLAGS_WARNINGS += -Wtautological-overlap-compare
		CXXFLAGS_WARNINGS += -Wtautological-type-limit-compare
		CXXFLAGS_WARNINGS += -Wtautological-unsigned-char-zero-compare
		CXXFLAGS_WARNINGS += -Wtautological-unsigned-enum-zero-compare
		CXXFLAGS_WARNINGS += -Wtautological-unsigned-zero-compare
		CXXFLAGS_WARNINGS += -Wtautological-value-range-compare
		CXXFLAGS_WARNINGS += -Wthread-safety
		CXXFLAGS_WARNINGS += -Wtype-limits
		CXXFLAGS_WARNINGS += -Wunaligned-access
		CXXFLAGS_WARNINGS += -Wundef
		CXXFLAGS_WARNINGS += -Wundef-prefix
		CXXFLAGS_WARNINGS += -Wundefined-func-template
		CXXFLAGS_WARNINGS += -Wundefined-internal-type
		CXXFLAGS_WARNINGS += -Wundefined-reinterpret-cast
		CXXFLAGS_WARNINGS += -Wunguarded-availability
		CXXFLAGS_WARNINGS += -Wuninitialized
		CXXFLAGS_WARNINGS += -Wunnamed-type-template-args
		CXXFLAGS_WARNINGS += -Wunneeded-internal-declaration
		CXXFLAGS_WARNINGS += -Wunneeded-member-function
		# # unknown: CXXFLAGS_WARNINGS += -Wuninitialized-const-reference
		CXXFLAGS_WARNINGS += -Wunreachable-code
		CXXFLAGS_WARNINGS += -Wunreachable-code-aggressive
		CXXFLAGS_WARNINGS += -Wunused
		# # unknown: CXXFLAGS_WARNINGS += -Wunsafe-buffer-usage
		# # unknown: CXXFLAGS_WARNINGS += -Wunsafe-buffer-usage-in-container
		CXXFLAGS_WARNINGS += -Wunsupported-dll-base-class-template
		CXXFLAGS_WARNINGS += -Wunused
		CXXFLAGS_WARNINGS += -Wunused-command-line-argument
		CXXFLAGS_WARNINGS += -Wused-but-marked-unused
		CXXFLAGS_WARNINGS += -Wuser-defined-literals
		CXXFLAGS_WARNINGS += -Wvariadic-macros
		CXXFLAGS_WARNINGS += -Wvector-conversion
		CXXFLAGS_WARNINGS += -Wvla
		# # disabled: CXXFLAGS_WARNINGS += -Wweak-vtables
		CXXFLAGS_WARNINGS += -Wwrite-strings
		CXXFLAGS_WARNINGS += -Wzero-as-null-pointer-constant
		CXXFLAGS_WARNINGS += -Wzero-length-array
		CXXFLAGS_WARNINGS += 
		
	endif
 
endif


# In any case, remove the following warnings 

# # DISABLED: CXXFLAGS_WARNINGS += -Wno-missing-braces

CXXFLAGS_WARNINGS += -Wunused-variable
CXXFLAGS_WARNINGS += -Wno-unused-parameter
CXXFLAGS_WARNINGS += -Wno-unknown-pragmas
CXXFLAGS_WARNINGS += -Wno-vla

CXXFLAGS_WARNINGS += -Wno-double-promotion	
CXXFLAGS_WARNINGS += -Wno-conversion 		# mostly unimportant stuff 
CXXFLAGS_WARNINGS += -Wno-sign-compare      # numerous comparisons between signed and unsigned indices 
CXXFLAGS_WARNINGS += -Wno-float-equal		# numerous comparisons of float to zero 
CXXFLAGS_WARNINGS += -Wno-sign-conversion	# too many signed indices 

ifeq ($(FLAG_CXX),GCC)

# # DISABLED: CXXFLAGS_WARNINGS += -Wno-shadow
# # DISABLED: CXXFLAGS_WARNINGS += -Wno-stack-usage

else ifeq ($(FLAG_CXX),CLANG)

CXXFLAGS_WARNINGS += -Wno-defaulted-function-deleted
CXXFLAGS_WARNINGS += -Wno-gnu-zero-variadic-macro-arguments
CXXFLAGS_WARNINGS += -Wno-vla-extension 
CXXFLAGS_WARNINGS += -Wno-shorten-64-to-32

endif














###############################################
#                                             #
#          Debugging information              #
#                                             #
###############################################

CXXFLAGS_DEBUG :=
ifeq ($(FLAG_NO_DEBUGINFO),yes)
else
CXXFLAGS_DEBUG += -g
endif


###############################################
#                                             #
#          Profiling instrumentation          #
#                                             #
###############################################

ifeq ($(FLAG_DO_PROFILE),yes)

	CXXFLAGS_PROF:= -pg -fno-omit-frame-pointer 

endif


###############################################
#                                             #
#          Sanitizer instrumentation          #
#                                             #
###############################################

# There are several incompatibilities between sanitizers 
# thread cannot be combined with address and leak
# Warning: memory causes considerable slowdown

# NOTE: The following construction of the strings ensures that there are no spaces added.

ifeq ($(FLAG_DO_USE_SANITIZER),yes)

	SANITIZERS :=
	
	ifeq ($(FLAG_CXX),GCC)

		SANITIZERS :=$(SANITIZERS)pointer-compare,pointer-subtract,
		SANITIZERS :=$(SANITIZERS)undefined,

	else ifeq ($(FLAG_CXX),CLANG)

		SANITIZERS :=$(SANITIZERS)pointer-compare,pointer-subtract,
		SANITIZERS :=$(SANITIZERS)undefined,

		SANITIZERS :=$(SANITIZERS)float-divide-by-zero,
		SANITIZERS :=$(SANITIZERS)unsigned-integer-overflow,
		SANITIZERS :=$(SANITIZERS)implicit-conversion,
		SANITIZERS :=$(SANITIZERS)nullability-arg,
		SANITIZERS :=$(SANITIZERS)nullability-assign,
		SANITIZERS :=$(SANITIZERS)nullability-return,

		#SANITIZERS :=$(SANITIZERS)memory,

	endif

	SANITIZERS :=$(SANITIZERS)address,leak,
	# SANITIZERS :=$(SANITIZERS)thread 

	# thread cannot be combined with address and leak 

	
	# SANITIZERS :=-ftrapv,

	# Comment out the following line to disable ALL built-in sanitizers 
	CXXFLAGS_SANI := -fsanitize=$(SANITIZERS) -ftrapv 

endif


###############################################
#                                             #
#    Use TCMalloc instead of std allocators   #
#                                             #
###############################################

ifeq ($(FLAG_USE_TCMALLOC),yes)
	CXXFLAGS_MALLOC=-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
else
	CXXFLAGS_MALLOC=
endif








 
 
###############################################
#                                             #
#          Static analysis                    #
#                                             #
###############################################

ifeq ($(FLAG_DO_STATICANALYSIS),yes)

	ifeq ($(FLAG_CXX),GCC)

		CXXFLAGS_STATICANALYSER := -fanalyzer -Wanalyzer-too-complex

	else ifeq ($(FLAG_CXX),CLANG)

		CXXFLAGS_STATICANALYSER := 

	endif

endif


###############################################
#                                             #
#        Format of diagnostic settings        #
#                                             #
###############################################

ifeq ($(FLAG_CXX),CLANG)

  CXXFLAGS_DIAGNOSISFORMAT := -fdiagnostics-show-template-tree

endif








###############################################
#                                             #
#     Include the overwrite file again,       #
#     this time the strings will be set       #
#                                             #
###############################################

# Use this file to overwrite the default settings above on a local machine
# At this point, also the compiler flags will be set 
-include $(projectdir)/OVERWRITE.COMPILE.mk



###############################################
#                                             #
#   CXXFLAGS - gather C++ compiler flags      #
#                                             #
###############################################

CXXFLAGS := 
CXXFLAGS += ${CXXFLAGS_LANG}
CXXFLAGS += ${CXXFLAGS_DIAGNOSISFORMAT}
CXXFLAGS += ${CXXFLAGS_WARNINGS}
CXXFLAGS += ${CXXFLAGS_STATICANALYSER}
CXXFLAGS += ${CXXFLAGS_DEBUG}
CXXFLAGS += $(CXXFLAGS_PROF)
CXXFLAGS += $(CXXFLAGS_SANI)
CXXFLAGS += ${CXXFLAGS_MALLOC}
CXXFLAGS += ${CXXFLAGS_OPTIMIZE}
CXXFLAGS += ${CXXFLAGS_CODEGEN}

CXXFLAGS := $(strip $(CXXFLAGS))


CXXFLAGS_EXECUTABLE:=
CXXFLAGS_EXECUTABLE+=$(CXXFLAGS)

ifeq ($(FLAG_DO_OPTIMIZE),yes)
	ifeq ($(FLAG_CXX),ICC)
		CXXFLAGS_EXECUTABLE+=-fwhole-program -fuse-linker-plugin
	else 
		ifeq ($(FLAG_CXX),GCC)
			CXXFLAGS_EXECUTABLE+=-fwhole-program -fuse-linker-plugin
		endif
	endif
else
endif



###############################################
#                                             #
#   CPPFLAGS - gather preprocessor flags      #
#                                             #
###############################################

CPPFLAGS := 


CURRENT_COMMIT_ID := $(shell git rev-parse HEAD)

CPPFLAGS += -D GIT_COMMIT_ID=\"$(CURRENT_COMMIT_ID)\"


ifeq ($(FLAG_DISABLE_CHECK_MESHES),yes)
CPPFLAGS += -DDO_NOT_CHECK_MESHES
endif

ifeq ($(FLAG_USE_BACKTRACER),yes)
CPPFLAGS += -DUSE_BACKTRACER
endif

ifeq ($(FLAG_DISABLE_ASSERTIONS),yes)
CPPFLAGS += -DNDEBUG
endif

ifeq ($(FLAG_DISCARD_ASSERT_MESSAGES),yes)
CPPFLAGS += -DDISCARD_ASSERT_MESSAGES
endif

ifeq ($(FLAG_USE_ORIGINAL_ASSERT_MACRO),yes)
CPPFLAGS += -DUSE_ORIGINAL_ASSERT_MACRO
endif

ifeq ($(FLAG_USE_PRIMITIVE_LOGGING),yes)
CPPFLAGS += -DUSE_PRIMITIVE_LOGGING
endif

ifeq ($(FLAG_COLORED_OUTPUT),yes)
CPPFLAGS += -DUSE_COLORED_OUTPUT
endif



ifeq ($(FLAG_DISABLE_STDLIBDEBUG),yes)
else 
CPPFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_ASSERTIONS -D_GLIBCXX_SANITIZE_VECTOR
endif

ifeq ($(FLAG_DO_USE_EXTENDED_PRECISION),yes)
CPPFLAGS += -DEXTENDED_PRECISION
endif
ifeq ($(FLAG_DO_USE_SINGLE_PRECISION),yes)
CPPFLAGS += -DSINGLE_PRECISION
endif

ifeq ($(FLAG_DO_OPTIMIZE),yes)
CPPFLAGS += -DELIDE_HOT_FUNCTIONS
endif

CPPFLAGS := $(strip $(CPPFLAGS))


###############################################
#                                             #
#         LDLIBS - gather linker flags        #
#                                             #
###############################################

ifeq ($(FLAG_USE_TCMALLOC),yes)
	LDLIBS :=-l:libtcmalloc.so.4
else
	LDLIBS :=
endif

LDLIBS := $(strip $(LDLIBS))



###############################################
#                                             #
#         Choose the type of linking          #
#         - static                            #
#         - dynamic                           #
#         - unspecified                       #
#                                             #
###############################################


LINKINGTYPE:=unspecified

ifeq ($(FLAG_DO_OPTIMIZE),yes)
	LINKINGTYPE:=objectfile
else
ifeq ($(OS),Windows_NT)
	LINKINGTYPE:=static
else 
	LINKINGTYPE:=dynamic
endif
endif

# LINKINGTYPE:=static















###############################################
#                                             #
#             Show all parameters             #
#                                             #
###############################################

# print all the compilation flags set manually or automatically
.PHONY: parameters
parameters:
	$(info FLAG_CXX                       = $(FLAG_CXX) )
	$(info FLAG_DISABLE_ASSERTIONS        = $(FLAG_DISABLE_ASSERTIONS) ) 
	$(info FLAG_USE_ORIGINAL_ASSERT_MACRO = $(FLAG_USE_ORIGINAL_ASSERT_MACRO) )
	$(info FLAG_DISABLE_CHECK_MESHES      = $(FLAG_DISABLE_CHECK_MESHES) ) 
	$(info FLAG_USE_PRIMITIVE_LOGGING     = $(FLAG_USE_PRIMITIVE_LOGGING) )
	$(info FLAG_DISABLE_STDLIBDEBUG       = $(FLAG_DISABLE_STDLIBDEBUG) ) 
	$(info FLAG_NO_EXCEPTIONS             = $(FLAG_NO_EXCEPTIONS) ) 
	$(info FLAG_DO_USE_EXTENDED_PRECISION = $(FLAG_DO_USE_EXTENDED_PRECISION) ) 
	$(info FLAG_DO_USE_SINGLE_PRECISION   = $(FLAG_DO_USE_SINGLE_PRECISION) ) 
	$(info FLAG_EXCESSIVE_WARNINGS        = $(FLAG_EXCESSIVE_WARNINGS) ) 
	$(info FLAG_DO_USE_SANITIZER          = $(FLAG_DO_USE_SANITIZER) ) 
	$(info FLAG_DO_STATICANALYSIS         = $(FLAG_DO_STATICANALYSIS) ) 
	$(info FLAG_DO_OPTIMIZE               = $(FLAG_DO_OPTIMIZE) ) 
	$(info FLAG_ENABLE_OPENMP             = $(FLAG_ENABLE_OPENMP) ) 
	$(info FLAG_USE_TCMALLOC              = $(FLAG_USE_TCMALLOC) ) 
	$(info FLAG_NO_DEBUGINFO              = $(FLAG_NO_DEBUGINFO) ) 
	$(info FLAG_USE_BACKTRACER            = $(FLAG_USE_BACKTRACER) ) 
	$(info FLAG_DO_STRIP                  = $(FLAG_DO_STRIP) ) 
	$(info FLAG_DO_PROFILE                = $(FLAG_DO_PROFILE) ) 
	@true
	$(info ) 
	@true
	$(info $(PATH))
	# echo $(PATH)
	# echo $$(PATH)
	# echo %PATH%
	$(info MAKE                     = $(MAKE) )
	$(info CXX                      = $(CXX) )
	$(info CXXFLAGS_LANG            = $(CXXFLAGS_LANG) ) 
	$(info CXXFLAGS_OPTIMIZE        = $(CXXFLAGS_OPTIMIZE) )
	$(info CXXFLAGS_CODEGEN         = $(CXXFLAGS_CODEGEN) ) 
	$(info CXXFLAGS_WARNINGS        = $(CXXFLAGS_WARNINGS) )
	$(info CXXFLAGS_DEBUG           = $(CXXFLAGS_DEBUG) ) 
	$(info CXXFLAGS_PROF            = $(CXXFLAGS_PROF) ) 
	$(info SANITIZERS               = $(SANITIZERS) ) 
	$(info CXXFLAGS_SANI            = $(CXXFLAGS_SANI) ) 
	$(info CXXFLAGS_MALLOC          = $(CXXFLAGS_MALLOC) ) 
	$(info CXXFLAGS_STATICANALYSER  = $(CXXFLAGS_STATICANALYSER) ) 
	$(info CXXFLAGS_DIAGNOSISFORMAT = $(CXXFLAGS_DIAGNOSISFORMAT) ) 	
	$(info CXXFLAGS_EXECUTABLE      = $(CXXFLAGS_EXECUTABLE) ) 
	$(info CPPFLAGS                 = $(CPPFLAGS) )
	$(info LDFLAGS                  = $(LDFLAGS) ) 
	$(info LDLIBS                   = $(LDLIBS) ) 
	$(info LINKINGTYPE              = $(LINKINGTYPE) ) 
	@true
