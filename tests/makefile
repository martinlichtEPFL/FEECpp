SHELL = /bin/sh

.PHONY: default 
default: build

.PHONY: build 
# build is set up below 

projectdir :=..

include ../common.compile.mk 

components:=
components+=basic
components+=utility
components+=combinatorics
components+=operators
components+=dense
components+=sparse
components+=solver
components+=mesh
components+=vtk
components+=fem
components+=solverfem
components+=sullivan2D
components+=sullivan3D
components+=whitney2D
components+=whitney3D

include tests.affices.mk

testsdir:=.
projectdir:=../
pathvar:=$(CURDIR)/../

contextdir:=./basic
context:=basic
include tests.rules.mk
#include tests.run.mk
build: basic.tests 

contextdir:=./utility
context:=utility
include tests.rules.mk
#include tests.run.mk
build: utility.tests

contextdir:=./combinatorics
context:=combinatorics
include tests.rules.mk
#include tests.run.mk
build: combinatorics.tests

contextdir:=./operators
context:=operators
include tests.rules.mk
#include tests.run.mk
build: operators.tests

contextdir:=./dense
context:=dense
include tests.rules.mk
#include tests.run.mk
build: dense.tests

contextdir:=./sparse
context:=sparse
include tests.rules.mk
#include tests.run.mk
build: sparse.tests

contextdir:=./solver
context:=solver
include tests.rules.mk
#include tests.run.mk
build: solver.tests

contextdir:=./mesh
context:=mesh
include tests.rules.mk
#include tests.run.mk
build: mesh.tests

contextdir:=./vtk
context:=vtk
include tests.rules.mk
#include tests.run.mk
build: vtk.tests

contextdir:=./fem
context:=fem
include tests.rules.mk
#include tests.run.mk
build: fem.tests

contextdir:=./solverfem
context:=solverfem
include tests.rules.mk
#include tests.run.mk
build: solverfem.tests

contextdir:=./sullivan2D
context:=sullivan2D
include tests.rules.mk
#include tests.run.mk
build: sullivan2D.tests

contextdir:=./sullivan3D
context:=sullivan3D
include tests.rules.mk
#include tests.run.mk
build: sullivan3D.tests

contextdir:=./whitney2D
context:=whitney2D
include tests.rules.mk
#include tests.run.mk
build: whitney2D.tests

contextdir:=./whitney3D
context:=whitney3D
include tests.rules.mk
#include tests.run.mk
build: whitney3D.tests





.PHONY: clean
clean:
	@echo "Finished cleaning test directory." 


.PHONY: dependclean
dependclean:
	@echo "Finished cleaning test directory dependency information files." 


.PHONY: outputclean
outputclean:
	@rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk ./*.svg ./*/*.svg ./*/*/*.svg
	@echo "Finished cleaning test directory output files."






	




