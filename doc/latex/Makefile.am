PROJECT = snalys
ROBOOPTS = C SORT
SRCDIR = @srcdir@

ALLFILES = $(notdir $(wildcard @top_srcdir@/src/*.f)) \
	$(notdir $(wildcard @top_srcdir@/src/*.inc))
# LATEX
#
LATEXDOCS=$(ALLFILES:=.tex)
LATEXAUX=$(ALLFILES:=.aux)
LATEXXREFS=$(LATEXDOCS:.tex=.tex.xref)
LATEXXREFSFILE=$(PROJECT)_tex.xrefs

%.tex.xref : @top_srcdir@/src/%
	-ln -f $< .
	@rdocpath@ $(notdir $<) $(notdir $(@:.xref=)) $(ROBOOPTS) INTERNAL GENXREF $@

# master index file, currently works only for html and latex documentation.
# Note that you can define the title of the document.
$(PROJECT)_mi.tex: $(LATEXXREFSFILE)
	@rdocpath@ $< $@ INDEX LATEX TITLE "$(PROJECT) API Reference"

# create xrefs file (file with the names of all .xref files).
$(LATEXXREFSFILE) : $(LATEXXREFS)
	/bin/ls  $(LATEXXREFS) > $@

# Rule to create latex documentation from a source file.
# We do not include source items, and generate laxtex documents
# than can be included in a master document.
%.tex : @top_srcdir@/src/%
	-ln -f $< .
	@rdocpath@ $(notdir $<) $(notdir ${@}) LATEX $(ROBOOPTS) NOSOURCE SINGLEDOC XREF $(LATEXXREFSFILE)

tex: $(LATEXXREFSFILE) $(LATEXDOCS) $(PROJECT)_mi.tex

# Use the latex programs to generate a .dvi from the master index file
# for our prokect. View this .dvi file with xdvi
ps :  tex
	latex $(PROJECT)_mi
	makeindex $(PROJECT)_mi
	latex $(PROJECT)_mi
	latex $(PROJECT)_mi
	dvips $(PROJECT)_mi.dvi -o


all : ps

clean :
	-rm *.tex *.log *.aux *.idx *.dvi *.xref *.xrefs \
	*.toc *.f *.inc *.ilg *.ind *.ps
distclean :
	-rm Makefile

# EOF
