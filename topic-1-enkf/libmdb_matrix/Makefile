CC = gcc
LIBTOOL = libtool
AR = ar

PWD := $(shell pwd)
UNAME := $(shell uname)
MACHINE := $(shell uname -m)

LIBS = -lm -lgsl
LDFLAGS = -L.

ifeq ($(UNAME),Linux)
   DEFINES = -DLINUX -DOPENBLAS
   LIBS += -llapacke
   ifeq ($(MACHINE),x86_64)
      DEFINES += -DLONG_PTR
       LIBS += -lopenblas
   else ifeq ($(MACHINE),i686)
      LDFLAGS += -L/usr/lib/sse2 -L/usr/lib/sse2/atlas
      LIBS += -latlas -lf77blas -lcblas
   else ifeq ($(MACHINE),ppc)
      LDFLAGS += -L/usr/lib/altivec
      LIBS += -llapack_atlas -lcblas \
              -Wl,-rpath,/usr/lib/altivec
   else
      $(error Unknown machine type $(MACHINE))
   endif
else ifeq ($(UNAME),Darwin)
   DEFINES = -DOSX -DLONG_PTR -DVECLIB
   LIBS += -lc -framework Accelerate
   INCLUDE_DIR += -I/usr/local/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers
else
   $(error $(UNAME) is not supported!)
endif

ELEM_TEST = double

MILD_WARNINGS = -Wall -Wstrict-prototypes
# The strict warning setup below is used by the FreeBSD people
STRICT_WARNINGS = -W -Wall -ansi -pedantic -Wbad-function-cast -Wcast-align \
                  -Wcast-qual -Wchar-subscripts -Winline \
                  -Wmissing-prototypes -Wnested-externs -Wpointer-arith \
                  -Wredundant-decls -Wshadow -Wstrict-prototypes -Wwrite-strings

CPPFLAGS = $(MILD_WARNINGS) $(INCLUDE_DIR)
CFLAGS = $(CPPFLAGS) -O3 -g
CFLAGS_PIC := $(CFLAGS) -fPIC

#DEFINES += $(DEFINES) -DNDEBUG

SO_LIBS_S := $(LIBS) -lfftw3f
SO_LIBS_D := $(LIBS) -lfftw3

ifeq ($(ELEM_TEST),float)
   DEFINES += -DFLOAT_ELEM_TEST
   LIBS += -lmdb_matrix_s -lfftw3f
else
ifeq ($(ELEM_TEST),double)
   DEFINES += -DDOUBLE_ELEM_TEST
   LIBS += -lmdb_matrix_d -lfftw3
else
endif
endif


#####################################################################


LIB_TARGET = libmdb_matrix_s.so libmdb_matrix_d.so

STATIC_LIB_TARGET = libmdb_matrix_s.a libmdb_matrix_d.a

TEST_TARGET = full_r_test full_c_test sparse_rcs_test sparse_rcs_mvm_test \
	      sparse_ccs_test toeplitz_test filter_test psf_test psf_2d_test \
              filter_new_test fftwe_test lapack_test llist_test \
	      psf_3d_test vector_test zfull_r_test sparse_rcs_mmm_test \
              lt_c_test toeplitz_new_test lt_r_test sparse_lil_test \
	      counter_test util_test

TEST_TARGET_OBJ = $(addsuffix .o, $(TEST_TARGET))

# The order matters because the headers for each of these files gets
# combined into one by build_header.py
LIB_BASENAME = elem util blas fftwe vector full_r sparse_rcs sparse_ccs full_c \
	       counter full_ut_c full_lt_c full_lt_r diag toeplitz \
               filter filter_new psf psf_2d psf_3d stop_watch \
               zfull_r zfull_c lapack sparse_coo toeplitz_new llist sparse_lil


LIB_OBJ_S = $(addsuffix _s.o, $(LIB_BASENAME))
LIB_OBJ_D = $(addsuffix _d.o, $(LIB_BASENAME))
LIB_H = boolean.h
LIB_H += $(addsuffix .h, $(LIB_BASENAME))


all: $(LIB_TARGET) $(TEST_TARGET) TAGS

ifeq ($(UNAME),Linux)
   libmdb_matrix_s.so: mdb_matrix_s.h
   libmdb_matrix_s.so: $(LIB_OBJ_S)
	$(CC) -shared -Wl,-soname,$(PWD)/$@ -o $@ $(LIB_OBJ_S) $(SO_LIBS_S) $(LDFLAGS)
else ifeq ($(UNAME),Darwin)
   libmdb_matrix_s.so: mdb_matrix_s.h
   libmdb_matrix_s.so: $(LIB_OBJ_S)
	$(LIBTOOL) -dynamic -o $@ $(LIB_OBJ_S) -install_name $(PWD)/$@ $(SO_LIBS_S) $(LDFLAGS)
endif

# The order of the header files appear matters!
mdb_matrix_s.h: $(LIB_H)
	./build_header.py -a "#define FLOAT_ELEM" -o $@ $^


ifeq ($(UNAME),Linux)
   libmdb_matrix_d.so: mdb_matrix_d.h
   libmdb_matrix_d.so: $(LIB_OBJ_D)
	$(CC) -shared -Wl,-soname,$(PWD)/$@ -o $@ $(LIB_OBJ_D) $(SO_LIBS_D) $(LDFLAGS)

else ifeq ($(UNAME),Darwin)
   libmdb_matrix_d.so: mdb_matrix_d.h
   libmdb_matrix_d.so: $(LIB_OBJ_D)
	$(LIBTOOL) -dynamic -o $@ $(LIB_OBJ_D) -install_name $(PWD)/$@ $(SO_LIBS_D) $(LDFLAGS)
endif

# The order the header files appear matters!
mdb_matrix_d.h: $(LIB_H)
	./build_header.py -a "#define DOUBLE_ELEM" -o $@ $^


static: $(STATIC_LIB_TARGET)

libmdb_matrix_s.a: mdb_matrix_s.h
libmdb_matrix_s.a: $(LIB_OBJ_S)
libmdb_matrix_s.a:
	$(AR) rcs $@ $(LIB_OBJ_S)

libmdb_matrix_d.a: mdb_matrix_d.h
libmdb_matrix_d.a: $(LIB_OBJ_D)
libmdb_matrix_d.a:
	$(AR) rcs $@ $(LIB_OBJ_D)


full_r_test: full_r_test.o

full_c_test: full_c_test.o

zfull_r_test: zfull_r_test.o

sparse_rcs_test: sparse_rcs_test.o

sparse_rcs_mvm_test: sparse_rcs_mvm_test.o

sparse_ccs_test: sparse_ccs_test.o

toeplitz_test: toeplitz_test.o

vector_test: vector_test.o

filter_new_test: filter_new_test.o

filter_test: filter_test.o

psf_test: psf_test.o

psf_2d_test: psf_2d_test.o

psf_3d_test: psf_3d_test.o

sparse_rcs_mmm_test: sparse_rcs_mmm_test.o

lt_c_test: lt_c_test.o

lt_r_test: lt_r_test.o

toeplitz_new_test: toeplitz_new_test.o

fftwe_test: fftwe_test.o

lapack_test: lapack_test.o

llist_test: llist_test.o

sparse_lil_test: sparse_lil_test.o

counter_test: counter_test.o

util_test: util_test.o


$(TEST_TARGET): mdb_matrix_d.h mdb_matrix_s.h

$(TEST_TARGET_OBJ): mdb_matrix_d.h mdb_matrix_s.h
$(TEST_TARGET): libmdb_matrix_s.so libmdb_matrix_d.so


TAGS: *.[ch]
	etags *.[ch]

%_s.o : %.c
	$(CC) -c $(CFLAGS_PIC) $(DEFINES) -DFLOAT_ELEM $< -o $@


%_d.o : %.c
	$(CC) -c $(CFLAGS_PIC) $(DEFINES) -DDOUBLE_ELEM $< -o $@

%.o : %.c
	$(CC) -c $(CFLAGS) $(DEFINES) $< -o $@

% : %.o
	$(CC) $(LDFLAGS) -o $@ $< $(LIBS)


################################################################################
# Adapted from Section 4.14 Generating Prerequisites Automatically of
# the reference manual for GNU make

DEPEND_S = $(patsubst %.c,%_s.d,$(wildcard *.c))
DEPEND_D = $(patsubst %.c,%_d.d,$(wildcard *.c))

%_s.d: %.c
	@set -e; rm -f $@; \
        $(CC) -MM -MG -MT $(patsubst %.d,%.o,$@) $(CPPFLAGS) $(DEFINES) -DFLOAT_ELEM $< > $@.$$$$; \
        sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
        rm -f $@.$$$$

%_d.d: %.c
	@set -e; rm -f $@; \
        $(CC) -MM -MG -MT $(patsubst %.d,%.o,$@) $(CPPFLAGS) $(DEFINES) -DDOUBLE_ELEM $< > $@.$$$$; \
        sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
        rm -f $@.$$$$

include $(DEPEND_S) $(DEPEND_D)

################################################################################

clean:
	rm -f *.o *.d $(LIB_TARGET) mdb_matrix_s.h mdb_matrix_d.h $(TEST_TARGET)
