#
# This makefile targets the installation of architecture-independent data.
#

include ../Makefile.inc

install: spelling.txt superatom.txt
	$(INSTALL_DIR) $(DESTDIR)$(datadir)
	$(INSTALL_DATA) $? $(DESTDIR)$(datadir)

uninstall:
	$(RM) -f $(DESTDIR)$(datadir)/spelling.txt
	$(RM) -f $(DESTDIR)$(datadir)/superatom.txt
	
../Makefile.inc: ../Makefile.inc.in ../config.status
	cd .. && ./config.status
