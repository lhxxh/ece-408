CC = gcc

UNAME := $(shell uname)
MACHINE := $(shell uname -m)

LIBS = -lm -lgsl -lconfuse

PATH_MDB_MATRIX ?= ../libmdb_matrix#$(HOME)/src/libmdb_matrix


ifeq ($(UNAME),Linux)
   #DEFINES = -DATLAS -DLINUX
   DEFINES = -DLINUX -DOPENBLAS
   LDFLAGS = #-Wl,-rpath=$(PATH_MDB_MATRIX)

   ifeq ($(MACHINE),x86_64)
      INCLUDE_DIR = -I$(PATH_MDB_MATRIX) #-I/usr/local/atlas/include
      LDFLAGS += -L$(PATH_MDB_MATRIX) #-L/usr/local/atlas/lib
      #LIBS += -latlas -lptf77blas -lptcblas
      LIBS += -lopenblas
   else ifeq ($(MACHINE),i686)
      LDFLAGS += -L/usr/lib/sse2 -L/usr/lib/sse2/atlas
      LIBS += -latlas -lf77blas -lcblas
   else ifeq ($(MACHINE),ppc)
      LDFLAGS += -L$(PATH_MDB_MATRIX) -L/usr/lib/altivec -Wl,-rpath,/usr/lib/altivec
      LIBS += -latlas -lf77blas -lcblas
   else
      $(error Unknown machine type $(MACHINE))
   endif
else ifeq ($(UNAME),Darwin)
   DEFINES = -DOSX -DLONG_PTR -DVECLIB
   LIBS += -lc -framework Accelerate
   INCLUDE_DIR += -I/usr/local/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers
   INCLUDE_DIR += -I$(PATH_MDB_MATRIX)
   LDFLAGS = -L$(PATH_MDB_MATRIX)
else
   $(error $(UNAME) not supported!)
endif

ELEM = double

MILD_WARNINGS = -Wall -Wstrict-prototypes
# The strict warning setup below is used by the FreeBSD people
STRICT_WARNINGS = -W -Wall -ansi -pedantic -Wbad-function-cast -Wcast-align \
                  -Wcast-qual -Wchar-subscripts -Winline \
                  -Wmissing-prototypes -Wnested-externs -Wpointer-arith \
                  -Wredundant-decls -Wshadow -Wstrict-prototypes -Wwrite-strings
CPPFLAGS = $(INCLUDE_DIR) #$(MILD_WARNINGS) 
CFLAGS = $(CPPFLAGS) -O3 -g

ifeq ($(ELEM),float)
   DEFINES += -DLENKF_FLOAT_ELEM
   LIBS += -lfftw3f -lmdb_matrix_s
else
ifeq ($(ELEM),double)
   DEFINES += -DLENKF_DOUBLE_ELEM
   LIBS += -lfftw3 -lmdb_matrix_d
else
endif
endif

# Uncomment the following to enamble profiling via gprof
#CFLAGS += -pg
#LDFLAGS += -pg

#DEFINES += $(DEFINES) -DNDEBUG

#DEFINES += -DHASH_FUNCTION=HASH_FNV -DHASH_EMIT_KEYS=3


################################################################################


BIN = compute_P_HT compute_P_HT_cuda#lenkf ensemble_test randn_test blas_test randn_conv_test_new convmtx

all: $(BIN) #TAGS

ensemble_test: ensemble.o ensemble_test.o randn.o util.o

randn_test: randn_test.o randn.o ensemble.o util.o

randn_conv_test_new: randn_conv_test_new.o randn.o ensemble.o

lenkf: lenkf.o lenkf_main.o ensemble.o randn.o util.o \
       lenkf_config.o arg_bundle.o edot_table.o

blas_test: blas_test.o randn.o util.o

convmtx: convmtx.o

compute_P_HT: compute_P_HT.o randn.o ensemble.o lenkf.o \
              arg_bundle.o edot_table.o util.o

compute_P_HT_cuda: compute_P_HT_cuda.cu
	nvcc $< -o $@

TAGS: *.[ch]
	etags -a ${PATH_MDB_MATRIX}/TAGS *.[ch]


indent:
	$(INDENT) $(INDENT_FLAGS) *.c *.h

%.o : %.c
	$(CC) -c $(CFLAGS) $(DEFINES) $< -o $@

% : %.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LIBS)


################################################################################
# Adapted from Section 4.14 Generating Prerequisites Automatically of
# the reference manual for GNU make

DEPEND = $(patsubst %.c,%.d,$(wildcard *.c))

%.d: %.c
	@set -e; rm -f $@; \
        $(CC) -MM $(CPPFLAGS) $(DEFINES) $< > $@.$$$$; \
        sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
        rm -f $@.$$$$

include $(DEPEND)

################################################################################


clean:
	rm -f *.o *.d $(BIN) TAGS
