#!/usr/bin/make

default: build

projectdir   :=..
include $(projectdir)/common.compile.mk

moddir :=.
module:=$(notdir $(CURDIR))

include $(projectdir)/common.module.mk

build: $(module).build