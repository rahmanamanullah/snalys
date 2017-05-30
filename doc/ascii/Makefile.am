ROBOOPTS = C SORT

SRCDIR = @srcdir@

ALLFILES = $(notdir $(wildcard @top_srcdir@/src/*.f)) \
	$(notdir $(wildcard @top_srcdir@/src/*.inc))
ASCIIDOCS=$(ALLFILES:=.txt)
%.txt : @top_srcdir@/src/%
	@rdocpath@ $< $@ ASCII 
ascii: $(ASCIIDOCS)
all : ascii
clean :
	-rm *.txt
distclean : clean
	-rm Makefile
