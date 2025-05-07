TARGETS = \
CahnHilliard2D 

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


CahnHilliard2D: CahnHilliard2D.o
	$(CLINKER) -o $@ $< $(PETIGA_LIB)
	$(RM) -f $<

OPTS=-nox -malloc_debug -malloc_dump

include $(PETIGA_DIR)/lib/petiga/conf/test
