bin_PROGRAMS = ct-cbn h-cbn
dist_doc_DATA = CHANGELOG README

# ct-cbn
ct_cbn_SOURCES = src/ct-cbn.c src/queue.c src/ct-cbn.h src/queue.h
ct_cbn_CFLAGS = $(GSL_CFLAGS) $(OPENMP_CFLAGS)
ct_cbn_LDADD = $(GSL_LIBS) $(OPENMP_CFLAGS)

# h-cbn
h_cbn_SOURCES = src/h-cbn.c src/queue.c src/ct-cbn.h src/queue.h
h_cbn_CFLAGS = $(GSL_CFLAGS) $(OPENMP_CFLAGS)
h_cbn_LDADD = $(GSL_LIBS) $(OPENMP_CFLAGS)

# example files
dist_pkgdata_DATA = \
	examples/BC.label \
	examples/BC.pat \
	examples/BC.poset \
	examples/CRC.label \
	examples/CRC.pat \
	examples/CRC.poset \
	examples/HIV.gps \
	examples/RCC.label \
	examples/RCC.pat \
	examples/RCC.poset \
	examples/hiv.pat \
	examples/hiv.poset \
	examples/hiv.prf \
	examples/prostate.pat \
	examples/prostate.poset \
	examples/prostate.prf
