SUBDIRS = src doc
dist_doc_DATA = README

.PHONY: test valtest
check:
	cd $(top_builddir)/tests &&	make check

test:
	cd $(top_builddir)/tests && $(MAKE) $(AM_MAKEFLAGS) test

valtest:
	cd $(top_builddir)/src && $(MAKE) $(AM_MAKEFLAGS) valtest
