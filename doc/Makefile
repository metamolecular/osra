#
# This makefile targets the compilation and installation of poject documentation.
#
# We provide compiled documentation on purpose, to reduce the number of dependencies (see linux/debian/control).
#

include ../Makefile.inc

HTML_DOC_DIR=html
NATURALDOCS_DIR=nd

all:

ifdef NATURALDOCS
all: $(HTML_DOC_DIR)/index.html

$(HTML_DOC_DIR)/index.html:
	mkdir -p $(HTML_DOC_DIR) $(NATURALDOCS_DIR)
	$(NATURALDOCS) -i ../src -o HTML $(HTML_DOC_DIR) -p $(NATURALDOCS_DIR)
endif

ifdef XSLTPROC
all: $(NAME).1 $(NAME).html

$(NAME).1: manual.sgml
#	This command has been taken from http://www.debian.org/doc/manuals/maint-guide/ch-dother.en.html#s-manpagexml:
	$(XSLTPROC) --nonet \
		--param make.year.ranges 1 \
		--param make.single.year.ranges 1 \
		--param man.charmap.use.subset 0 \
		http://docbook.sourceforge.net/release/xsl/current/manpages/docbook.xsl \
		$?

$(NAME).html: manual.sgml
	$(XSLTPROC) --nonet \
		--param make.year.ranges 1 \
		--param make.single.year.ranges 1 \
		--param man.charmap.use.subset 0 \
		http://docbook.sourceforge.net/release/xsl/current/html/docbook.xsl \
		$? > $@

install: $(NAME).1
	$(INSTALL_DIR) $(DESTDIR)$(mandir)/man1
	$(INSTALL_DATA) $? $(DESTDIR)$(mandir)/man1
endif

uninstall:
	$(RM) -f $(DESTDIR)$(mandir)/man1/$(NAME).1

clean:
	$(RM) -Rf $(HTML_DOC_DIR) $(NAME).1 $(NAME).html

distclean:
	$(RM) -Rf $(NATURALDOCS_DIR) manual.sgml

../Makefile.inc: ../Makefile.inc.in ../config.status
	cd .. && ./config.status
