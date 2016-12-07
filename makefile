IFE3DMAIN_DIR = .
BIN_DIR = $(IFE3DMAIN_DIR)

#.SUFFIXES= .f .c .o
PROG_VIBRA = ../bin/vibra
PROG_VIBRA_P = ../bin/vibra_presion

OBJECTS = lecmalla.o signop.o signof.o\
write_interfaz_fs.o write_interfaz_fp.o write_interfaz_sp.o\
read_interfaz_fs.o read_interfaz_fp.o read_interfaz_sp.o\
ifindkey.o face2elem.o volus.o puntar.o puntars.o\
putitap.o intersec.o intersec2.o ordnod.o\
putitas.o intersecs.o ordnods.o putitasp.o\
matrpl.o matep1.o matep2.o ensp.o\
coefp.o coefs.o densip.o densis.o densif.o velso.o\
matrso.o mates1.o mates2.o enss.o\
matrfl.o matef1.o matef2.o ensf.o ccoipmos.o\
coipmo.o ccoipmo.o ccoipmo2.o coismo.o ccoismo.o coifmo.o\
coipmos.o coismos.o blomap.o blomas.o blomaf.o robin.o\
reaeltyp.o getsol.o cnum.o trasimf.o vec_alloc.o\
domain_div.o vtkcolor.o blomaf_presion.o getsol_presion.o\
matrfl_presion.o matef1_presion.o matef2_presion.o\
depre_presion.o grad_presion.o vtkfl_presion.o\
vtkfl_case_presion.o\
matrflpl_presion.o ensflpl_presion.o\
matrflso_presion.o ensflso_presion.o\
puntar_presion.o intersec_presion.o ensflpl_presion2.o\
neighbors.o metis_perm.o metis_nodend_f2c.o\
read_grafo.o write_grafo.o create_grafo.o
METISOBJS=ndmetis_f2c.o smbfactor.o cmdline_ndmetis.o io.o\
computefillin_f2c.o
EOBJECTS = caljee.o trasim.o conel4.o\
cal_componentesp.o cal_componentesf.o solver4x4.o\
cramer3x3.o cal_jlf.o cal_uhf.o elem2ref.o searstr.o depre.o\
contfl.o diver.o valor.o vtkso.o vtkpl.o vtkfl.o vwmerrfl.o\
vtkso_case.o vtkpl_case.o vtkfl_case.o vwmpl.o vwmso.o\
solp_conf.o solp_noconf.o co_rayleigh.o co_rayleigh2.o\
desplazamientofluido_vtk.o desplazamientofluido_vtk_case.o

CFLAGS      = -Wall $(USER_INCLUDE)
#FPPFLAGS    = $(USER_INCLUDE)
FFLAGS       = -Wno-tabs #$(USER_INCLUDE)
VERSION         = 2.3
BASE_DIR        = $(PAR_IFE_DIR)
#LOCDIR          = $(IFE3DMAIN_DIR)
#HDIR            = $(BASE_DIR)/src/par_gpfep/headers
#ODIR            = $(GPMAIN_DIR)/obj/$(P_ARCH)
#SDIR            = $(BASE_DIR)/src/par_gpfep
#GPFEP_DIR       = $(BASE_DIR)/src/gpfep
#GPLIB_DIR       = $(BASE_DIR)/lib/$(P_ARCH)
#PAR_GPLIB_DIR   = $(GPLIB_DIR)
#USER_INCLUDE    = -I$(HDIR)
#PROF            = 
#HEADERS         = 
#SOURCES         = 
METISGKLIB       = $(METISDIR)/GKlib
METISINCL        = $(METISDIR)/include
METISLIB         = $(METISDIR)/libmetis
METISPROG        = $(METISDIR)/programs

all: $(METISOBJS) vibra vibra_presion 

vibra: $(OBJECTS) $(EOBJECTS) vibra.o 
		-$(FLINKER) -o $(PROG_VIBRA) vibra.o $(OBJECTS) $(EOBJECTS) $(METISOBJS)\
		$(PETSC_FORTRAN_LIB) $(PETSC_LIB) -lz ${SLEPC_LIB} ${METIS_LIB}

vibra_presion: $(OBJECTS) $(EOBJECTS) vibra_presion.o
		-$(FLINKER) -o $(PROG_VIBRA_P) vibra_presion.o $(OBJECTS) $(EOBJECTS) $(METISOBJS)\
		$(PETSC_FORTRAN_LIB) $(PETSC_LIB) -lz ${SLEPC_LIB} ${METIS_LIB}

vibra.o: vibra.F cfgreading.f

vibra_presion.o: vibra_presion.F cfgreading.f

ndmetis_f2c.o: ndmetis_f2c.c 
		$(PETSC_COMPILE_SINGLE) -I$(METISGKLIB) -I$(METISINCL) \
		-I$(METISLIB) -I$(METISPROG) $< 

cmdline_ndmetis.o: cmdline_ndmetis.c 
		$(PETSC_COMPILE_SINGLE) -I$(METISGKLIB) -I$(METISINCL) \
		-I$(METISLIB) -I$(METISPROG) $< 

smbfactor.o: smbfactor.c 
		$(PETSC_COMPILE_SINGLE) -I$(METISGKLIB) -I$(METISINCL) \
		-I$(METISLIB) -I$(METISPROG) $< 

io.o: io.c 
		$(PETSC_COMPILE_SINGLE) -I$(METISGKLIB) -I$(METISINCL) \
		-I$(METISLIB) -I$(METISPROG) $< 

computefillin_f2c.o: computefillin_f2c.c 
		$(PETSC_COMPILE_SINGLE) -I$(METISGKLIB) -I$(METISINCL) \
		-I$(METISLIB) -I$(METISPROG) $< 

clean_:
		rm -f *.o

wipe: clean
		rm -f $(PROGRAM) *% *~

include $(SLEPC_COMMON)
