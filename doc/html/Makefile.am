PROJECT = snalys
SRCDIR = @srcdir@

ALLFILES = $(notdir $(wildcard @top_srcdir@/src/*.f)) \
	$(notdir $(wildcard @top_srcdir@/src/*.inc))

# The various documentation files, derived from the source files.
# HTML	
#
HTMLDOCS=$(ALLFILES:=.html)
HTMLXREFS=$(HTMLDOCS:.html=.html.xref)
HTMLXREFSFILE=$(PROJECT)_html.xrefs

# Rule to create an .xref file from a source file for the various formats.
%.html.xref : @top_srcdir@/src/%
	@rdocpath@ $< $(notdir $(@:.xref=)) $(ROBOOPTS) INTERNAL GENXREF $@

# master index file, currently works only for html and latex documentation.
# Note that you can define the title of the document.
$(PROJECT)_mi.html: $(HTMLXREFSFILE) 
	@rdocpath@ $< $@ INDEX HTML TITLE "$(PROJECT) Master Index"

# create xrefs file (file with the names of all .xref files).
$(HTMLXREFSFILE) : $(HTMLXREFS)
	/bin/ls $(HTMLXREFS) > $@

# Rule to create html documentation from a source file.
%.html : @top_srcdir@/src/%
	@rdocpath@ $< $(notdir ${@}) HTML $(ROBOOPTS) XREF $(HTMLXREFSFILE)

# Create the xref files for the various formats.
xhtml: $(HTMLXREFSFILE) 

# Create the documentation files for the various formats.
html: $(HTMLXREFSFILE) $(HTMLDOCS) $(PROJECT)_mi.html

all : html

clean :
	-rm *.html *.xref *.xrefs

distclean : clean
	-rm Makefile

# EOF
