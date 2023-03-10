# Warning: dependences on .h files are not given. If you change
# these, type "make clean; make"

OBJS=check2xsf.o check_read.o xsf_write.o molecule_fix.o cube_write.o \
     cell_read.o periodic_table.o basis.o super.o xplor_write.o \
     pdb_write.o cell_write.o dx_write.o vasp_write.o pdb_read.o \
     xyz_write.o cml_write.o fdf_write.o gpfa.o rotate.o chdiff_read.o \
     esp_read.o interpolate.o line_write.o shelx_read.o shelx_write.o \
     fbin_write.o ident_sym.o primitive.o ksym.o f15.o cif_write.o \
     cube_read.o cif_read.o c2x2spg.o sort_atoms.o py_write.o \
     xsf_read.o vasp_read.o denfmt_write.o denfmt_read.o parse.o \
     dipole.o fort34_read.o crystal_read.o abinit_read.o abinit_write.o \
     abinit_in_read.o fdf_read.o qe_write.o qe_read.o parity.o potential.o \
     qe_rho_read.o qe_xml_read.o dict.o charge.o qe_psi_read.o print_cell.o \
     bands_write.o geom_write.o xv_read.o rho_read.o tube.o ccp4_write.o \
     xc.o band_process.o gcoeff_write.o wavecar_write.o bxsf_write.o \
     bands_read.o xyz_read.o gcoeff_read.o elk_write.o elk_read.o \
     geom_read.o npy_read.o npy_write.o file_read.o data_combine.o \
     dx_read.o


# Recommended: -DSPGLIB

# For SPGLIB (1.8.3 or greater)

DEFS=-DSPGLIB
CPPFLAGS=-I.
# Add "-fopenmp" to next line if the linker complains about various
# missing threads
LDFLAGS=-L.
LIBS=-lsymspg

# To compile without SPGLIB, set all of the above to be empty

#DEFS=
#CPPFLAGS=
#LDFLAGS=
#LIBS=

# Linux or MacOS / x86_64 / gcc
CFLAGS=-Wall -Wno-unused-result -O -g $(DEFS)

# Linux / IA32 / gcc
#CFLAGS=-Wall -Wno-unused-result -O -D_FILE_OFFSET_BITS=64 $(DEFS) -g
# Tru64
#CFLAGS=-O -std1
# Solaris
#CFLAGS=-O -xarch=native64

CC=gcc

c2x: c2xsf.h $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o  c2x $(OBJS) -lm $(LIBS)

check2xsf.o: c2xsf.h

clean:
	-rm c2x $(OBJS)

