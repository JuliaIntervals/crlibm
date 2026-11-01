## Process this file with automake to produce Makefile.in


lib_LIBRARIES = libscs.a
INCLUDES=-I$(top_srcdir)

if COMPILER_ICC
AM_CFLAGS =  -std=c99 -mp -Qoption,cpp,--extended_float_types -IPF_fp_speculationsafe
else
if USE_SSE2
AM_CFLAGS = -std=c99 -Wall -O3 -msse2 -mfpmath=sse -ffp-contract=off
else
AM_CFLAGS = -std=c99 -Wall -O3 -ffp-contract=off
endif
endif


libscs_a_SOURCES = scs.h scs_private.h scs_private.c \
	double2scs.c print_scs.c \
	division_scs.c\
	addition_scs.c multiplication_scs.c\
	scs2double.c zero_scs.c \
	scs2mpfr.c 


include_HEADERS = scs.h ../crlibm_config.h

