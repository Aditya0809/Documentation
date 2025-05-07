TARGETS = \
ch2d_implicit \
ch2d_explicit \ 

ALL: $(TARGETS)
clean::
	-@$(RM) $(TARGETS)

CFLAGS    = #-g3 -Wall -Wextra -Wno-unused-parameter #-Wconversion
FFLAGS    = #-g3 -Wall -Wextra -fcheck=all
CPPFLAGS  =
FPPFLAGS  =
LOCDIR    = demo/
EXAMPLESC =
EXAMPLESF =
MANSEC    = IGA

topdir := $(shell cd .. && pwd)
PETIGA_DIR ?= $(topdir)
include $(PETIGA_DIR)/lib/petiga/conf/variables
include $(PETIGA_DIR)/lib/petiga/conf/rules


ch2d_implicit: ch2d_implicit.o
	$(CLINKER) -o $@ $< $(PETIGA_LIB)
	$(RM) -f $<
ch2d_explicit: ch2d_explicit.o
	$(CLINKER) -o $@ $< $(PETIGA_LIB)
	$(RM) -f $<

OPTS=-nox -malloc_debug -malloc_dump

include $(PETIGA_DIR)/lib/petiga/conf/test
