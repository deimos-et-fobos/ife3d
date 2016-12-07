      PROGRAM PAR_IFE3D
   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
     
C PETSC include
!#include "finclude/petsc.h"
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscis.h"
#include "finclude/petscdraw.h" 
#include "finclude/petscviewer.h"
C SLEPC includes
!#include "finclude/slepc.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"
#include "finclude/slepcst.h"
C

      PARAMETER(NDELP=120000,NDVERP=100000,NDARP=200000,NDGLP=
     &      5*NDVERP+NDARP,NDELS=10000000,NDVERS=3000000,NDGLS=3*NDVERS,
     &      NDELF=5000000,NDVERF=5000000,NDNODF=10000000,NDGLF=NDNODF,
     &      NDCACOF=1500000,ndvercos=200000,ndvercof=83000,
     &      NDUMITAS=5000000,NDUMITAP=5000000,
     &      NDDOF=NDGLP+NDGLS+NDGLF+NDCACOF,MAXSD=30,MAXREF=100,
     &      ndproc=100)      
      

      DIMENSION MMP(3,NDELP),NNP(6,NDELP),NREP(NDELP),NRAP(3,NDELP),
     &      NRVP(3,NDELP),NSDP(NDELP),ZP(2,NDVERP),ISDP(MAXSD),
     &      youngp(MAXSD),poisp(MAXSD),corrkp(MAXSD),DENSP(MAXSD),
     &      IREFCP(MAXREF),NGP(3,NDELP),IFLAGP(NDARP),SMP(NDELP),
     &      THIN(NDELP)
      
      DIMENSION MMS(4,NDELS),NRVS(4,NDELS),NSDS(NDELS),ZS(3,NDVERS),
     &      ISDS(MAXSD),youngs(MAXSD),poiss(MAXSD),DENSS(MAXSD),
     &      NRCS(4,NDELS),IREFCS(MAXREF),VOLS(NDELS),iflag_vers(NDVERS)

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VOLF,VFSUM,
     &      DIV,PF,P_CONF,uf,iver
      INTEGER, DIMENSION(:), ALLOCATABLE :: NSDF,IFLAGF,vec_perm,
     &      vec_iperm,ietf,jetf,jeef
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: UFV,
     &      COMPONENTES_F,zf,nuf
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MMF,NRCF,NNF,NGF,
     &      FACE2EL,IETF2
      DIMENSION ISDF(MAXSD),DENSF(MAXSD),VELSQ(MAXSD),
     &      IREFCF(MAXREF),IREFRF(MAXREF),RK(MAXREF),RM(MAXREF)
!      DIMENSION MMF(4,NDELF),NRCF(4,NDELF),NNF(4,NDELF),
!     &      NSDF(NDELF),ZF(3,NDVERF),VOLF(NDELF),vfsum(ndnodf),
!     &      ISDF(MAXSD),DENSF(MAXSD),VELSQ(MAXSD),
!     &      IREFCF(MAXREF),IREFRF(MAXREF),RK(MAXREF),RM(MAXREF),
!     &      NGF(4,NDELF),IFLAGF(NDNODF),uf(ndglf)

!      INTEGER FACE2EL(5,NDNODF)
      INTEGER DIRP(NDELP),DOF,proj_pf,proj_sf,proj_sp,div_perm
      INTEGER ALLOCATESTATUS,DEALLOCATEStatus,CNUMP(NDVERP+NDARP)
      INTEGER color_p(ndverp),color_s(ndvers),color_f(ndverf)
!      DOUBLE PRECISION NUF(3,NDNODF)
      
      DIMENSION IVERCO(3,NDCACOF),ICREF(NDCACOF),XA(NDCACOF),
     &            ICACO(NDCACOF)

      DIMENSION ivercos(ndvercos),ielps(ndvercos),indlp(ndvercos)
      DIMENSION ivercof(ndvercof),ielpf(ndvercof),indlpf(ndvercof)

      DIMENSION UP(NDGLP),US(NDGLS),PRES(NDCACOF),PRESSP(NDVERCOS)
!      DIMENSION DIV(NDELF),P(NDVERF),UFV(3,NDVERF),IVER(NDVERF)
      DIMENSION UPV(3,NDVERP),BETA_L(2,NDVERP),ROT(2,NDARP+NDVERP)
      DIMENSION USV(3,NDVERS)

      INTEGER w_interfaz_fs,w_interfaz_fp,w_interfaz_sp
      INTEGER r_interfaz_fs,r_interfaz_fp,r_interfaz_sp
      INTEGER write_matrix,read_matrix
      integer r_mgrafo,w_mgrafo,r_ggrafo,w_ggrafo
      
      CHARACTER*80 STR,MALLAP,MALLAS,MALLAF,WorkDir,filegd,TITLE,
     &             graph_name
      CHARACTER*20 VARF
      character*2  Str2

!     DE LA RUTINA PUTITAS
      DIMENSION INODFS(NDUMITAS),ILCIFS(NDUMITAS),IELS(NDUMITAS),
     &            X1S(NDUMITAS),X2S(NDUMITAS),X3S(NDUMITAS),
     &            Y1S(NDUMITAS),Y2S(NDUMITAS),Y3S(NDUMITAS),
     &            Z1S(NDUMITAS),Z2S(NDUMITAS),Z3S(NDUMITAS),
     &            INDCS(NDUMITAS)

!     DE LA RUTINA PUTITAP
      DIMENSION INODFP(NDUMITAP),ILCIFP(NDUMITAP),IELP(NDUMITAP),
     &            X1P(NDUMITAP),X2P(NDUMITAP),X3P(NDUMITAP),
     &            Y1P(NDUMITAP),Y2P(NDUMITAP),Y3P(NDUMITAP),
     &            ZINTP(NDUMITAP),INDCP(NDUMITAP)

!      dimension ietf(ndverf+1),jetf(4*ndelf),jeef(4*ndelf)
      dimension iets(ndvers+1),jets(4*ndels),jees(4*ndels)
      dimension ietp(ndverp+ndarp+1),jetp(6*ndelp)
!      dimension ethaf(ndelf),ietf2(2,ndnodf)
!      dimension ethau(ndelf)
!      dimension componentes_f(4,ndelf)
      dimension componentes_pz(3,ndelp),componentes_prot(2,3,ndelp)
      dimension ethaf_aux(ndverf)
      dimension ndiv(3),cmin(3),cmax(3)

!*    NEIGHBORS, GRAFO & METIS_PERM
      integer, dimension(:), allocatable :: ia,ja,iat,jat,iat_v,jat_v,
     &  perm, iperm, ia_perm, ia_iperm  
      
! De la solucion conforme y no-conforme de la presión del fluido
!      dimension p_conf(ndverf),p_noconf(ndnodf),zf_noconf(3,ndnodf)

! Para cada parte del estimador en el fluido
!      double precision jlf(ndelf),uhf(ndelf)
      
      REAL*4 TI,TA(2)
      
      PetscLogDouble mem
      PetscViewer Mviewer,Kviewer
      PetscScalar EIG_R,EIG_I,IOVECSCAN(1)
      PetscScalar SHIFT,SHIFT2
      PetscOffset I_X
      PetscReal   TOL,EIG_ERROR,TOLERR,TOLRES
      PetscReal   cutoff
      PetscInt    NEV,EPSMAXIT,ITS,NCONV,ICONT
      PetscInt    MAXIT,its_bal
      PetscInt    r_min,r_max
      PetscInt    alloc_d(NDDOF)
      PetscInt    alloc_o(NDDOF)
      PetscInt    pc_nl,pc_fl
      MPI_Comm    COMM 
      Vec         IOVEC    ! VECTOR FOR INPUT/OUTPUT (R)
      Vec         B        ! RHSS
      Vec         X        ! UNKNOWNS VECTORS OF AX=B
      Mat         KR       ! STIFFNESS MATRICES
      Mat         M        ! MASS MATRICES
      EPS         eps      ! EIGENVALUE PROBLEM SOLVER!
      EPSBalance  eps_bal_type
      VecScatter  PAR2SEQ  ! SCATTER FROM PARALLEL TO SEQUENTIAL
      KSP         ksp,subksp(100)      ! KRYLOV SUBSPACE CONTEXTS
      PC          pc,subpc       ! PRECONDITIONER CONTEXTS
      ST          st       ! SPECTRAL TRANSFORMATION CONTEXT
      CHARACTER*20 PCTYPE  ! PRECONDITIONER TYPE (CHAR *)
      CHARACTER*20 KSPTYPE ! KRYLOV SUBSPACE METHOD TYPE (CHAR *)
      CHARACTER*20 epstype ! EPS METHOD TYPE (CHAR *)
      CHARACTER*20 sttype  ! ST METHOD TYPE (CHAR *)
      EPSProblemType epsproblemtype ! EPS PROBLEM TYPE (CHAR *)
      CHARACTER*40  EMSG04
      CHARACTER*80  EMSG08
      CHARACTER*75  ERRMSG
      CHARACTER*120 EMSG12
      
      common/coelasp/youngp,poisp,corrkp,densp,ndsdp,isdp
      COMMON/coelass/youngs,poiss,denss,ndsds,isds
      COMMON/COEFFL/DENSF,VELSQ,NDSDF,ISDF
      common/p_fluido/varf
      EXTERNAL COEFP,COEFS,DENSIF,VELSO
      
!*
!*********************************************************************
!*    INSERT HERE ANY DEFINITIONS OF ADDITIONAL VARIABLES OR ARRAYS
!*********************************************************************
!*
!*
!*********************************************************************
!*    INITIALIZING
!******************************************************************** 
!*
      CALL SlepcInitialize (PETSC_NULL_CHARACTER,IERR)
      CALL MPI_Comm_size (PETSC_COMM_WORLD,ISIZE,IERR)
      CALL MPI_Comm_rank (PETSC_COMM_WORLD,IRANK,IERR)
      
      TTIME0 = MPI_Wtime()
      WRITE(*,*) ISIZE, IRANK
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      
      CALL PetscPrintf (PETSC_COMM_WORLD,
     &      CHAR(10)//' STARTING IFE3D'//CHAR(10)//CHAR(10), IERR)
      CALL PetscPrintf (PETSC_COMM_WORLD,
     &      ' -------------------------------'//CHAR(10), IERR)

!*********************************************************************
!     Allocamiento de memoria
!*********************************************************************
      ALLOCATE ( volf(NDELF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** volf ***"
      ALLOCATE ( vfsum(ndnodf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** vfsum ***"
      ALLOCATE ( nuf(3,ndnodf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** nuf ***"
      ALLOCATE ( ZF(3,NDVERF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** zf ***"
      ALLOCATE ( div(NDELF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** div ***"
      ALLOCATE ( pf(NDVERF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** pf ***"
      ALLOCATE ( ufv(3,NDVERF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** ufv ***"
      ALLOCATE ( uf(ndglf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** uf ***"
      ALLOCATE ( iver(NDVERF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** iver ***"
      ALLOCATE ( componentes_f(4,ndelf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** componentes_f ***"
      ALLOCATE ( p_conf(ndverf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** p_conf ***"
      ALLOCATE ( mmf(4,NDELF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** mmf ***"
      ALLOCATE ( nrcf(4,NDELF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** nrcf ***"
      ALLOCATE ( nnf(4,NDELF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** nnf ***"
      ALLOCATE ( nsdf(NDELF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** nsdf ***"
      ALLOCATE ( ngf(4,NDELF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** ngf ***"
      ALLOCATE ( iflagf(NDNODF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** iflagf ***"
      ALLOCATE ( face2el(5,NDNODF), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** face2el ***"
      ALLOCATE ( ietf2(2,ndnodf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** ietf2 ***"
      ALLOCATE ( ietf(ndverf+1), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** ietf ***"
      ALLOCATE ( jetf(4*ndelf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** jetf ***"
      ALLOCATE ( jeef(4*ndelf), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** jeef ***"
!*********************************************************************


!*
!******     ENTRADA DE DATOS - SOLO EL PROCESADOR PRINCIPAL
!*
      IF(IRANK.EQ.0) THEN     ! HASTA CONSTRUCION DE KR y M
      
         IF (ndproc.LT.isize) THEN
            PRINT*,' TOO MANY PROCESSES ---> NDPROC = ',ndproc
            STOP
         ENDIF

         TI=DTIME(TA)
         WRITE(*,*) "Reading of general data: Started"
         filegd = 'ife3d.cfg'
         include 'cfgreading.f'
         WRITE(*,*) "Reading of general data: Complete"
!*
!******  LECTURA DE MALLAS.
!*
         WRITE(*,*) "Reading meshes..."
         CALL LECMALLA(MALLAP,MALLAS,MALLAF,INDPLA,INDSOL,INDFLU,
     &            DIRP,SMP,THIN,NELP,NARP,NVERP,MMP,NNP,NRAP,NREP,
     &            NRVP,ZP,NSDP,NELS,NVERS,MMS,NRVS,NRCS,ZS,NSDS,
     &            NELF,NNODF,NVERF,NNF,MMF,NRCF,ZF,NSDF)

         NGLP=5*NVERP+NARP
         NGLS=3*NVERS
         NGLF=NNODF
         DOF=NGLP+NGLS+NGLF
         varf='desplazamientos'
         
         TI=DTIME(TA)
         WRITE(*,*) "---> " ,TI, " SEC."

!*
!******     Chequeo de dimensiones.
!*
         IF (NDVERP.LT.NVERP.OR.NDELP.LT.NELP.OR.NDGLP.LT.NGLP
     &                                       .OR.NDARP.LT.NARP) THEN
            PRINT*,' UNDERSIZED'
            PRINT*,' NDVERP= ',NDVERP,' NVERP= ',NVERP
            PRINT*,' NDARP= ',NDARP,' NARP= ',NARP
            PRINT*,' NDELP= ',NDELP,' NELP= ',NELP
            PRINT*,' NDGLP= ',NDGLP,' NGLP= ',NGLP
            STOP
         ENDIF

         IF (NDVERS.LT.NVERS.OR.NDELS.LT.NELS.OR.NDGLS.LT.NGLS) THEN
            PRINT*,' UNDERSIZED'
            PRINT*,' NDVERS= ',NDVERS,' NVERS= ',NVERS
            PRINT*,' NDELS= ',NDELS,' NELS= ',NELS
            PRINT*,' NDGLS= ',NDGLS,' NGLS= ',NGLS
            STOP
         ENDIF

         IF (NDVERF.LT.NVERF.OR.NDELF.LT.NELF.OR.NDGLF.LT.NGLF) THEN
            PRINT*,' UNDERSIZED'
            PRINT*,' NDVERF= ',NDVERF,' NVERF= ',NVERF
            PRINT*,' NDELF= ',NDELF,' NELF= ',NELF
            PRINT*,' NDGLF= ',NDGLF,' NGLF= ',NGLF
            STOP
         ENDIF

!*
!******     DETERMINACION DE LOS VECTORES NORMALES EN LAS CARAS 
!*               DEL FLUIDO. DEVUELVE EN NGF(L,K) 1 O -1
!*
         IF (INDPLA.EQ.1) THEN
            TI=DTIME(TA)
            WRITE(*,*) "PLATE RT sign determination"
            CALL SIGNOP(NNP,NELP,NARP+NVERP,IFLAGP,NGP)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
            
            ! CONSTRUCCION DEL PUNTERO CNUMP.
            TI=DTIME(TA)
            WRITE(*,*) "PLATE matrix pointers"
            CALL CNUM(NELP,NVERP,NARP,NNP,CNUMP)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
*            WRITE(12,*)(CNUMP(I),I=1,NVERP+NARP+1)
         ENDIF
         IF (INDFLU.EQ.1) THEN 
            TI=DTIME(TA)
            WRITE(*,*) "FLUID RT sign determination"
            CALL SIGNOF(NNF,NELF,NNODF,IFLAGF,NGF)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF

!*
!******  DETERMINACION DE LOS ELEMENTOS (T+) Y (T-) DE LAS CARAS
!*          Y LOS VECTORES NORMALES (T+) --> (T-) EN EL FLUIDO
!*
         IF (INDFLU.EQ.1) THEN
            TI=DTIME(TA)
            WRITE(*,*) "FLUID face/element connectivity"
            CALL FACE2ELEM(ZF,MMF,NNF,FACE2EL,NGF,NELF,NNODF,NUF,VOLF)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF

!*
!******  DETERMINACION DEL VOLUMEN DE LOS ELEMENTOS DEL SOLIDO 3D.
!*
         IF (INDSOL.EQ.1) THEN
            TI=DTIME(TA)
            WRITE(*,*) "SOLID elements volumes calculation"
            CALL VOLUS(ZS,MMS,NELS,NVERS,VOLS)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF

!*
!******     DETERMINACION DE LAS CARAS DE CONTACTO DEL FLUIDO.
!*
         IF (INDPLA*INDFLU.EQ.1 .OR. INDSOL*INDFLU.EQ.1) THEN
           WRITE(*,*) "PUNTAR (FLUID)"
           CALL PUNTAR(ZF,MMF,NNF,NELF,NRCF,NRCOF,IREFCF,
     &                 NCACOF,ncacofp,ICACO,IVERCO,ICREF,XA,
     &                 nrcop,nrcos,irefcp,irefcs)
           nvercof=0
           ncacofs=ncacof-ncacofp
           if(ncacofs.eq.0) then
             write(*,*)"---> No SOLID/FLUID interface detected"
           else
             write(*,*)"---> Potential SOLID/FLUID interface detected"
           endif
           if(ncacofp.eq.0) then
             write(*,*)"---> No PLATE/FLUID interface detected"
           else
             write(*,*)"---> Potential PLATE/FLUID interface detected"
           endif
           if(ncacof.gt.0) then
             write(*,*)"---> ",'NDCACOF =',NDCACOF,' NCACOF =',NCACOF
           endif

           TI=DTIME(TA)
           WRITE(*,*) "---> " ,TI, " SEC."
         ELSE
           NCACOF=0
           NCACOFP=0
           NCACOFS=0
           nvercof=0
         ENDIF

         !* CHEQUEO DE DIMENSIONES.
         IF (NDCACOF.LT.NCACOF) THEN
            PRINT*,' UNDERSIZED'
            PRINT*,' NDCACOF= ',NDCACOF,' NCACOF= ',NCACOF
            STOP
         ENDIF

!*
!******        ARMADO DE LAS CARITAS DE CONTACTO (SOLIDO/FLUIDO).
!*
         IF (INDSOL*INDFLU.EQ.1 .and. ncacofs.gt.0) THEN
          WRITE(*,*) "CARITAS (SOLID/FLUID)"
          if(r_interfaz_fs.eq.0)then 
            WRITE(*,*) "---> Computing SOLID/FLUID interface"
            CALL PUTITAS(NELS,NRCOS,NRCS,IREFCS,NCACOF,MMS,ZS,ZF,
     &               IVERCO,NUMITAS,INODFS,ICACO,ILCIFS,ICREF,IELS,
     &               INDCS,X1S,X2S,X3S,Y1S,Y2S,Y3S,Z1S,Z2S,Z3S,proj_sf)

            if(w_interfaz_fs.eq.1)then
              WRITE(*,*) "---> Writing interface_fs.dat"
              call write_interfaz_fs(numitas,inodfs,ilcifs,
     &                  iels,indcs,x1s,x2s,x3s,y1s,y2s,y3s,z1s,z2s,z3s)
            endif
          else
            WRITE(*,*) "---> Reading interface_fs.dat"
            call read_interfaz_fs(numitas,inodfs,ilcifs,
     &               iels,indcs,x1s,x2s,x3s,y1s,y2s,y3s,z1s,z2s,z3s)
          endif

          if(numitas.eq.0) then
            write(*,*)"---> No SOLID/FLUID interface detected"
          else
            write(*,*)"---> SOLID/FLUID interface detected"
            WRITE(*,*)"---> ",'NDUMITAS =',NDUMITAS,' NUMITAS =',NUMITAS
          endif

          TI=DTIME(TA)
          WRITE(*,*) "---> " ,TI, " SEC."

          !* CHEQUEO DE DIMENSIONES.
          IF (NDUMITAS.LT.NUMITAS) THEN  
            PRINT*,' UNDERSIZED'
            PRINT*,' NDUMITAS= ',NDUMITAS,' NUMITAS= ',NUMITAS
            STOP
          ENDIF

         else
            numitas = 0   
         ENDIF    !(INDSOL.AND.INDFLU)

!*
!******        ARMADO DE LAS CARITAS DE CONTACTO (PLACA/FLUIDO).
!*
         IF (INDPLA*INDFLU.EQ.1 .and. ncacofp.gt.0) THEN
          WRITE(*,*) "CARITAS (PLATE/FLUID)"
          if(r_interfaz_fp.eq.0)then
            WRITE(*,*) "---> Computing PLATE/FLUID interface"
            CALL PUTITAP(NELP,NRCOP,NRAP,NREP,IREFCP,NCACOF,MMP,ZP,ZF,
     &            IVERCO,NUMITAP,INODFP,ICACO,ILCIFP,ICREF,IELP,INDCP,
     &            X1P,X2P,X3P,Y1P,Y2P,Y3P,ZINTP,THIN,SMP,DIRP,proj_pf) 
!
            if(w_interfaz_fp.eq.1)then
              WRITE(*,*) "---> Writing interface_fp.dat"
              call write_interfaz_fp(numitap,inodfp,ilcifp,ielp,indcp,
     &                                  x1p,x2p,x3p,y1p,y2p,y3p,zintp)
            endif
          else
            WRITE(*,*) "---> Reading interface_fp.dat"
            call read_interfaz_fp(numitap,inodfp,ilcifp,ielp,indcp,
     &                               x1p,x2p,x3p,y1p,y2p,y3p,zintp)
          endif
!
          if(numitap.eq.0) then
            write(*,*)"---> No PLATE/FLUID interface detected"
          else
            write(*,*)"---> PLATE/FLUID interface detected"
            WRITE(*,*)"---> ",'NDUMITAP =',NDUMITAP,' NUMITAP =',NUMITAP
          endif
!
          TI=DTIME(TA)
          WRITE(*,*) "---> " ,TI, " SEC."
         else
            numitap = 0   
         ENDIF    
         !* CHEQUEO DE DIMENSIONES.
         IF (NDUMITAP.LT.NUMITAP) THEN
           PRINT*,' UNDERSIZED'
           PRINT*,' NDUMITAP= ',NDUMITAP,' NUMITAP= ',NUMITAP
           STOP
         ENDIF

!*
!******     DETERMINACION DE LAS CARAS DE CONTACTO DEL SOLIDO.
!*
         IF (INDPLA*INDSOL.EQ.1) THEN
           WRITE(*,*) "PUNTAR (SOLID)"
           if(r_interfaz_sp.eq.0)then
             WRITE(*,*) "---> Computing SOLID/PLATE interface"
             CALL PUNTARS(ZS,MMS,nvers,NELS,nrvs,NRCS,NRCOS,IREFCS,
     &                    nvercos,iflag_vers,zp,mmp,nelp,nrap,nrcop,
     &                    irefcp,thin,smp,dirp,ivercos,ielps,indlp)
             if(w_interfaz_sp.eq.1)then
               WRITE(*,*) "---> Writing interface_sp.dat"
               call write_interfaz_sp(nvercos,ivercos,ielps,indlp)
             endif
           else
             WRITE(*,*) "---> Reading interface_sp.dat"
             call read_interfaz_sp(nvercos,ivercos,ielps,indlp)
           endif
!
           if(nvercos.eq.0) then
             write(*,*)"---> No SOLID/PLATE interface detected"
           else
             write(*,*)"---> SOLID/PLATE interface detected"
             write(*,*)"---> ",'NDVERCOS =',NDVERCOS,
     &                         ' NVERCOS =',NVERCOS
           endif
!
           TI=DTIME(TA)
           WRITE(*,*) "---> " ,TI, " SEC."
         ELSE
           nvercos=0
         ENDIF
         !* CHEQUEO DE DIMENSIONES.
         IF (ndvercos.LT.nvercos) THEN
            PRINT*,' UNDERSIZED'
            PRINT*,' NDVERCOS= ',NDVERCOS,' NVERCOS= ',NVERCOS
            STOP
         ENDIF
!
         DOF=DOF+NCACOF+3*NVERCOS
         IF (DOF.GT.NDDOF) THEN
            WRITE(*,*) 'ERROR: DOF = ',DOF,'---> NDDOF = ',NDDOF
            STOP
         ENDIF

!*
!******    PERMUTACION DE COLUMNAS
!*
         ALLOCATE ( vec_perm(dof), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** vec_perm ***"
         ALLOCATE ( vec_iperm(dof), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** vec_iperm ***"
         DO 1115 i=1,dof
           vec_perm(i)=i-1    ! start with 0
           vec_iperm(i)=i-1   ! start with 0
 1115    continue
!
!*       Permutacion con METIS por default
         if(metis.eq.1) then
           WRITE(*,*) "COLUMN PERMUTATION WITH METIS"
           n_ia = nelp + nels + nelf + nvercos + numitap + numitas 
           n_ja = 6*nelp + 4*nels + 4*nelf + 8*nvercos + 
     &            8*numitap + 6*numitas 
           n_iat = nverp + narp + nvers + nnodf + nvercos +
     &             ncacofp + ncacofs
           n_jatv = 1
           ALLOCATE ( ia(n_ia+1), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** ia ***"
           ALLOCATE ( ja(n_ja), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** ja ***"
           ALLOCATE ( iat(n_iat+1), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** iat ***"
           ALLOCATE ( jat(n_ja), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** jat ***"
           ALLOCATE ( iat_v(n_iat+1), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** iat_v ***"
           ALLOCATE ( perm(n_iat), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** perm ***"
           ALLOCATE ( iperm(n_iat), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** iperm ***"
           ALLOCATE ( ia_perm(n_iat+1), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** ia_perm ***"
           ALLOCATE ( ia_iperm(n_iat+1), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** ia_iperm ***"
!
           if(r_mgrafo.eq.1)then
             graph_name="grafo_perm.dat"
             write(*,*) '---> Reading graph from '//trim(graph_name)
             call read_grafo(n_iat,nedge_grafo,iat_v,jat_v,
     &                       graph_name,r_mgrafo,1) 
             n_jatv = iat_v(n_iat+1)-1
             ALLOCATE (jat_v(n_jatv),STAT = AllocateStatus)
             IF (AllocateStatus /= 0) STOP "*** jat_v ***"
             call read_grafo(n_iat,nedge_grafo,iat_v,jat_v,
     &                       graph_name,r_mgrafo,2) 
           endif
           if(r_mgrafo.eq.0)then
             call neighbors(ia,ja,n_ia,n_ja,nelp,nels,nelf,nvercos,
     &          numitap,numitas,nverp,narp,nvers,nnodf,ncacofp,
     &          ncacofs,nnp,mms,nnf,ielps,ivercos,ielp,ilcifp,inodfp,
     &          iels,ilcifs,inodfs,iat,jat,iat_v,nedge_grafo)
             if(allocated(jat_v))then
               deallocate(jat_v,STAT=DEALLOCATEStatus)
               IF(DEALLOCATEStatus /= 0) STOP "ERROR in deallocation."
             endif
             n_jatv = iat_v(n_iat+1)-1
             ALLOCATE (jat_v(n_jatv),STAT = AllocateStatus)
             IF (AllocateStatus /= 0) STOP "*** jat_v ***"
             call create_grafo(ia,ja,iat,jat,iat_v,jat_v,n_iat)
             if(w_mgrafo.eq.1)then
               graph_name="grafo_perm.dat"
               write(*,*) '---> Writing graph in '//trim(graph_name)
               call write_grafo(n_iat,nedge_grafo,iat_v,jat_v,
     &                          graph_name)
             endif
           endif
           call metis_perm(iat_v,n_iat,nedge_grafo,nverp,narp,nvers,
     &          nnodf,nvercos,ncacofp,ncacofs,cnump,ia_perm,ia_iperm,
     &          jat_v,perm,iperm,vec_perm,vec_iperm)
           DEALLOCATE(ia,ja,iat,jat,iat_v,jat_v,perm,iperm,ia_perm,
     &                ia_iperm,STAT = DEALLOCATEStatus)
           IF(DEALLOCATEStatus /= 0) STOP "ERROR in deallocation." 
!           WRITE(*,*) "---> Fill-In Analysis:"
!           call computefillin_f2c(n_iat,nedge_grafo,iat_v,jat_v,
!     &                            vec_perm,vec_iperm)
!
!           PERMUTACION UTILIZANDO EL GRAFO GLOBAL
!           call metis_perm_g(n_iat,nedge_grafo,iat_v,jat_v,
!     &                       vec_perm,vec_iperm)
!           call computefillin_f2c(n_iat,nedge_grafo,iat_v,jat_v,
!     &                            vec_perm,vec_iperm)
!           DEALLOCATE(ia,ja,iat,jat,STAT = DEALLOCATEStatus)
!           IF(DEALLOCATEStatus /= 0) STOP "ERROR in deallocation: "
!
         endif
!
!*       Permutacion con division del dominio
         if(div_perm.eq.1)then
           WRITE(*,*) "COLUMN PERMUTATION WITH DOMAIN DIVISION"
           call domain_div(cmin,cmax,ndiv,nglp,nelp,zp,mmp,nnp,cnump,
     &          ngls,nvers,zs,mms,nglf,nelf,zf,mmf,nnf,ncacof,ncacofp,
     &          icaco,iverco,nvercos,ivercos,vec_perm,vec_iperm,nverf,
     &          color_p,color_s,color_f,dirp,smp)
!           call vtkcolor(color_p,color_s,color_f,nverp,nvers,nverf)
!           WRITE(*,*) "---> Fill-In Analysis:"
!           call computefillin_f2c(n_iat,nedge_grafo,iat_v,jat_v,
!     &                            vec_perm,vec_iperm)
         endif 
         TI=DTIME(TA)
         WRITE(*,*) "---> " ,TI, " SEC."

!*
!***     GRAFO GLOBAL DEL PROBLEMA
!*
         WRITE(*,*) "GLOBAL GRAPH OF THE PROBLEM"
         n_ia = nelp + nels + nelf + nvercos + numitap + numitas 
         n_ja = 18*nelp + 12*nels + 4*nelf + 24*nvercos + 
     &          20*numitap + 14*numitas 
         n_iat = dof
         n_jatv = 1
         ALLOCATE ( ia(n_ia+1), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** ia ***"
         ALLOCATE ( ja(n_ja), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** ja ***"
         ALLOCATE ( iat(n_iat+1), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** iat ***"
         ALLOCATE ( jat(n_ja), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** jat ***"
         ALLOCATE ( iat_v(n_iat+1), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** iat_v ***"
         if(r_ggrafo.eq.1)then
           graph_name="grafo_global.dat"
           write(*,*) '---> Reading graph from '//trim(graph_name)
           call read_grafo(n_iat,nedge_grafo,iat_v,jat_v,
     &                     graph_name,r_ggrafo,1) 
           n_jatv = iat_v(n_iat+1)-1
           ALLOCATE (jat_v(n_jatv),STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** jat_v ***"
           call read_grafo(n_iat,nedge_grafo,iat_v,jat_v,
     &                     graph_name,r_ggrafo,2) 
         endif
         if(r_ggrafo.eq.0)then
           call neighbors_g(ia,ja,n_ia,n_ja,nelp,nels,nelf,nvercos,
     &          numitap,numitas,nverp,narp,cnump,nvers,nnodf,ncacofp,
     &          ncacofs,nnp,mms,nnf,ielps,ivercos,ielp,ilcifp,inodfp,
     &          iels,ilcifs,inodfs,iat,jat,iat_v,nedge_grafo)
           if(allocated(jat_v))then
             deallocate(jat_v,STAT=DEALLOCATEStatus)
             IF(DEALLOCATEStatus /= 0) STOP "ERROR in deallocation."
           endif
           n_jatv = iat_v(n_iat+1)-1
           ALLOCATE (jat_v(n_jatv),STAT = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** jat_v ***"
           call create_grafo(ia,ja,iat,jat,iat_v,jat_v,n_iat)
           if(w_ggrafo.eq.1)then
             graph_name="grafo_global.dat"
             write(*,*) '---> Writing graph in '//trim(graph_name)
             call write_grafo(n_iat,nedge_grafo,iat_v,jat_v,graph_name)
           endif
         endif
         DEALLOCATE(ia,ja,iat,jat,STAT = DEALLOCATEStatus)
         IF(DEALLOCATEStatus /= 0) STOP "ERROR in deallocation." 
         TI=DTIME(TA)
         WRITE(*,*) "---> " ,TI, " SEC."
!
         WRITE(*,*) "------------------------------"
         WRITE(*,*) "| D.O.F. Plate: ", nglp,"|"
         WRITE(*,*) "| D.O.F. Solid: ", ngls,"|"
         WRITE(*,*) "| D.O.F. Fluid: ", nglf,"|"
         WRITE(*,*) "| D.O.F. I.FP.: ", ncacofp,"|"
         WRITE(*,*) "| D.O.F. I.FS.: ", ncacofs,"|"
         WRITE(*,*) "| D.O.F. I.SP.: ", 3*nvercos,"|"
         WRITE(*,*) "------------------------------"
         WRITE(*,*) "| TOTAL D.O.F.: ", dof,"|"
         WRITE(*,*) "------------------------------"
! 
         CALL VecCreateSeq(PETSC_COMM_SELF,DOF,IOVEC,IERR)
      ENDIF       ! (IRANK.EQ.0)
!      
      CALL PetscMemoryGetCurrentUsage(mem, IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,'**************'//
     & '**************************************'//CHAR(10), IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      write(*,99)'Process ',IRANK,': memory used = ',mem/(10**6),' MB.'
 99   format(a8,I2,a16,f8.2,a4)

!*
!*********************************************************************
!*             CREACION DE LAS MATRICES KR Y M
!*********************************************************************
!*
      CALL MPI_Bcast(nglp,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ngls,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nglf,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ncacof,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nvercos,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nvercof,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(DOF,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(numitap,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(numitas,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nelp,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nels,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nelf,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nverp,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nvers,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nverf,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(narp,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nnp,6*nelp,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(mms,4*nels,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(mmf,4*nelf,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nnf,4*nelf,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(icaco,ncacof,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(iverco,3*ncacof,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(dirp,nelp,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ilcifp,numitap,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ielp,numitap,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ilcifs,numitas,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(iels,numitas,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ivercos,nvercos,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ielps,nvercos,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(indlp,nvercos,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ivercof,nvercof,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(ielpf,nvercof,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(indlpf,nvercof,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(NEV,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(MAXIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(INDFLU,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(INDSOL,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(INDPLA,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(SHIFT,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(TOLRES,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(TOLERR,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(n_jatv,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      if(irank.ne.0)then
        ALLOCATE ( vec_iperm(dof), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** vec_iperm ***"
        ALLOCATE (iat_v(dof+1), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** iat_v ***"
        ALLOCATE (jat_v(n_jatv),STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** jat_v ***"
      endif
      CALL MPI_Bcast(vec_iperm,dof,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(iat_v,dof+1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(jat_v,n_jatv,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,'**************'//
     & '**************************************'//CHAR(10), IERR)
      CALL VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,DOF,X,IERR)
      CALL VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,DOF,B,IERR)
      CALL VecScatterCreateToZero(X,PAR2SEQ,IOVEC,IERR)
!      
      CALL MatCreate(PETSC_COMM_WORLD,KR,IERR)
      CALL MatCreate(PETSC_COMM_WORLD,M,IERR)
      CALL MatSetType(KR,MATSBAIJ,IERR)
      CALL MatSetType(M,MATSBAIJ,IERR)
      CALL MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,DOF,DOF,IERR)
      CALL MatSetSizes(KR,PETSC_DECIDE,PETSC_DECIDE,DOF,DOF,IERR)
!
      if(read_matrix.eq.1 .and. write_matrix.eq.0) goto 6543
!      
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' ALLOCATING MEMORY FOR MATRICES M & KR'//CHAR(10), IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
!      
      TIME0 = MPI_Wtime()
      do 1111 i=1,dof
         alloc_d(i) = 0
         alloc_o(i) = 0
 1111 continue
!
      CALL VecGetOwnershipRange(X,r_min,r_max,ierr)
      CALL allocKR(dof,iat_v,jat_v,vec_iperm,r_min,r_max,
     &             alloc_d,alloc_o)
      CALL MatSeqSBAIJSetPreallocation(KR,1,1,alloc_d,IERR)
      CALL MatMPISBAIJSetPreallocation(KR,1,1,alloc_d,1,alloc_o,IERR)
      CALL MatSetOption(KR,MAT_SPD,PETSC_TRUE,IERR)
      CALL MatSetUp(KR,IERR)
!
      do 1112 i=1,dof
        alloc_d(i) = 0
        alloc_o(i) = 0
 1112 continue      
!
      CALL allocM(dof,iat_v,jat_v,vec_iperm,r_min,r_max,
     &            nglp,ngls,nglf,alloc_d,alloc_o)
      CALL MatSeqSBAIJSetPreallocation(M,1,1,alloc_d,IERR)
      CALL MatMPISBAIJSetPreallocation(M,1,1,alloc_d,1,alloc_o,IERR)
      CALL MatSetOption(M,MAT_SPD,PETSC_TRUE,IERR)
      CALL MatSetUp(M,IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
!
      CALL PetscMemoryGetCurrentUsage(mem, IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,'**************'//
     & '**************************************'//CHAR(10), IERR)
      write(*,99)'Process ',IRANK,': memory used = ',mem/(10**6),' MB.'
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,'**************'//
     & '**************************************'//CHAR(10), IERR)
      TIME1 = MPI_Wtime()
      WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC." 
      CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)

!*
!*********************************************************************
!*             CALCULO DE LAS MATRICES KR Y M
!*********************************************************************
!*
      IF(IRANK.EQ.0) THEN 
!*
!******  CALCULO DE LA MATRIZ DE RIGIDEZ Y DE MASA DE LA PLACA.
!*  
         IF (INDPLA.EQ.1) THEN
            WRITE(*,*) "PLATE matrices M & KR calculation"
            CALL MATRPL(ZP,MMP,NNP,NGP,NELP,CNUMP,NSDP,COEFP,
     &                  NUP,THIN,M,KR,vec_iperm)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF
!*
!******  CALCULO DE LA MATRIZ DE RIGIDEZ Y DE MASA DEL SOLIDO-3D.
!*  
         IF (INDSOL.EQ.1) THEN
            WRITE(*,*) "SOLID matrices M & KR calculation"
            CALL MATRSO(ZS,MMS,NELS,NGLP,NSDS,COEFS,VOLS,
     &                  M,KR,vec_iperm)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF
!*
!******     CALCULO DE LA MATRIZ DE RIGIDEZ Y DE MASA DEL FLUIDO.
!*  
         IF (INDFLU.EQ.1) THEN
           WRITE(*,*) "FLUID matrices M & KR calculation"
           CALL MATRFL(ZF,NNF,MMF,NGF,NELF,NGLP,NGLS,VOLF,
     &                 NSDF,DENSIF,VELSO,NRCF,NRRF,IREFRF,
     &                 RK,RM,M,KR,vec_iperm)
           TI=DTIME(TA)
           WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF
!*
!******                 CALCULO DE MATRICES DE CONTACTO.
!*  
         IF (INDPLA*INDFLU.EQ.1 .OR. INDSOL*INDFLU.EQ.1) THEN
           WRITE(*,*)"FLUID interaction matrix calculation"
           CALL COIFMO(ICACO,NCACOF,XA,NGLP,NGLS,DOF,KR,vec_iperm)
           TI=DTIME(TA)
           WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF
!
         IF (INDPLA*INDFLU.EQ.1) THEN
           WRITE(*,*)"FLUID/PLATE interaction matrix calculation"
           CALL COIPMO(NCACOF,NUMITAP,INODFP,ILCIFP,INDCP,
     &                 IELP,NELP,NUF,MMP,NNP,ZP,CNUMP,X1P,Y1P,
     &                 X2P,Y2P,X3P,Y3P,DIRP,SMP,DOF,KR,vec_iperm)
           TI=DTIME(TA)
           WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF
!
         IF (INDSOL*INDFLU.EQ.1) THEN
            WRITE(*,*)"FLUID/SOLID interaction matrix calculation"
            CALL COISMO(NCACOF,NUMITAS,INODFS,ILCIFS,IELS,
     &                  NELS,NUF,MMS,ZS,x1s,y1s,z1s,x2s,y2s,z2s,
     &                  x3s,y3s,z3s,NGLP,DOF,KR,vec_iperm)
            TI=DTIME(TA)
            WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF
!
         IF (INDPLA*INDSOL.EQ.1) THEN
           WRITE(*,*)"SOLID/PLATE interaction matrix calculation"
           CALL COIPMOS(nvercos,IVERCOS,ielps,INDLP,zs,mmp,nnp,zp,
     &                  CNUMP,DIRP,SMP,DOF,KR,ncacof,vec_iperm)
           CALL COISMOS(nvercos,IVERCOS,ielps,dirp,NGLP,DOF,KR,
     &                  ncacof,vec_iperm)
           TI=DTIME(TA)
           WRITE(*,*) "---> " ,TI, " SEC."
         ENDIF
!
      ENDIF       !(IRANK.EQ.0)

!*
!*********************************************************************
!*             ENSAMBLAJE DE LAS MATRICES KR Y M
!*********************************************************************
!*
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' ASSEMBLING MATRICES M & KR'//CHAR(10), IERR)
      TIME0 = MPI_Wtime()
      CALL MatAssemblyBegin( M, MAT_FINAL_ASSEMBLY,IERR)
      CALL MatAssemblyBegin(KR, MAT_FLUSH_ASSEMBLY,IERR)
      CALL MatAssemblyEnd( M, MAT_FINAL_ASSEMBLY,IERR)
      CALL MatAssemblyEnd(KR, MAT_FLUSH_ASSEMBLY,IERR)
      TIME1 = MPI_Wtime()
      WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC." 
      CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
!      CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,
!     &                                 "matM_presion.dat",Mviewer,IERR)
!      CALL PetscViewerSetFormat(Mviewer,PETSC_VIEWER_ASCII_MATLAB,IERR)
!      CALL MatView( M, Mviewer,IERR)

!*
!*********************************************************************
!*             CONDICIONES DE BORDE SOBRE LA MATRIZ KR
!*********************************************************************
!*
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' SETTING BOUNDARY CONDITIONS'//CHAR(10), IERR)
      TIME0 = MPI_Wtime()
      IF (IRANK.EQ.0) THEN
        CALL BLOMAP(KR,CNUMP,NNP,NRVP,NRAP,NELP,vec_iperm)
        CALL BLOMAS(KR,MMS,NRVS,NELS,NGLP,vec_iperm)
        CALL BLOMAF(KR,NNF,NRCF,NELF,NGLP,NGLS,vec_iperm)
      ENDIF
      TIME1 = MPI_Wtime()
      WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC." 
      CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
      
!*
!*********************************************************************
!*             REENSAMBLAJE DE LAS MATRICES KR Y M
!*********************************************************************
!*
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' RE-ASSEMBLING STIFFNESS MATRIX'//CHAR(10), IERR)
      TIME0 = MPI_Wtime()
      CALL MatAssemblyBegin(KR, MAT_FINAL_ASSEMBLY,IERR)     
      CALL MatAssemblyEnd(KR, MAT_FINAL_ASSEMBLY,IERR)
      TIME1 = MPI_Wtime()
      WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC." 
      CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)     
!      CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,
!     &                                 "matK_presion.dat",Kviewer,IERR)
!      CALL PetscViewerSetFormat(Kviewer,PETSC_VIEWER_ASCII_MATLAB,IERR)
!      CALL MatView(KR, Kviewer,IERR)

!*
!*********************************************************************
!*             ESCRITURA DE LAS MATRICES KR Y M
!*********************************************************************
!*
      if(write_matrix.eq.1)then
        CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' WRITING MASS MATRIX'//CHAR(10), IERR)
        TIME0 = MPI_Wtime()
!        call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'matM',
!     &                             FILE_MODE_WRITE,Mviewer,ierr)
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'matM',
     &                            Mviewer,ierr)
       CALL PetscViewerSetFormat(Mviewer,PETSC_VIEWER_ASCII_COMMON,IERR)
        call MatView( M,Mviewer,ierr)
        TIME1 = MPI_Wtime()
        WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC."
        CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
!
        CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' WRITING STIFFNESS MATRIX'//CHAR(10), IERR)
        TIME0 = MPI_Wtime()
!        call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'matK',
!     &                             FILE_MODE_WRITE,Kviewer,ierr)
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'matK',
     &                            Kviewer,ierr)
       CALL PetscViewerSetFormat(Kviewer,PETSC_VIEWER_ASCII_COMMON,IERR)
        call MatView(KR,Kviewer,ierr)
        TIME1 = MPI_Wtime()
        WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC."
        CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
      endif

!*
!*********************************************************************
!*             LECTURA DE LAS MATRICES KR Y M
!*********************************************************************
!*
 6543 if(read_matrix.eq.1 .and. write_matrix.eq.0)then
        CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' READING MASS MATRIX'//CHAR(10), IERR)
        TIME0 = MPI_Wtime()
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'matM',
     &                             FILE_MODE_READ,Mviewer,ierr)
        call MatLoad( M,Mviewer,ierr)
        TIME1 = MPI_Wtime()
        WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC."
        CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
!
        CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' READING STIFFNESS MATRIX'//CHAR(10), IERR)
        TIME0 = MPI_Wtime()
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'matK',
     &                             FILE_MODE_READ,Kviewer,ierr)
        TIME1 = MPI_Wtime()
        WRITE (EMSG04,*)"---> " ,TIME1 - TIME0, " SEC."
        CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
        call MatLoad(KR,Kviewer,ierr)
      endif
!      call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,
!     &              PETSC_NULL_CHARACTER,0,0,600,600,Kviewer,ierr)
!      call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,
!     &              PETSC_NULL_CHARACTER,600,0,600,600,Mviewer,ierr)
!      call MatView(KR,Kviewer,ierr)
!      call MatView(M,Mviewer,ierr)

!*
!*********************************************************************
!*                         SOLVER EPS
!*********************************************************************
!*
      CALL PetscMemoryGetCurrentUsage(mem, IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,'**************'//
     & '**************************************'//CHAR(10), IERR)
      write(*,99)'Process ',IRANK,': memory used = ',mem/(10**6),' MB.'
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,'**************'//
     & '**************************************'//CHAR(10), IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      TIME0 = MPI_Wtime()
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     CHAR(10)//' SOLVING....'//CHAR(10), IERR)
!     
      CALL EPSCreate(MPI_COMM_WORLD, eps, IERR)
      CALL EPSSetOperators(eps, KR, M, IERR)
      CALL EPSSetProblemType(eps, EPS_GHEP, IERR)
      SHIFT2=abs(SHIFT)*SHIFT
      CALL EPSSetTarget(eps,SHIFT,IERR)
      CALL EPSSetTolerances(eps,TOLERR,MAXIT,IERR)
      CALL EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL, IERR)
      CALL EPSSetDimensions(eps, NEV, PETSC_DECIDE, PETSC_DECIDE,IERR)
      CALL EPSGetST(eps,st,IERR)
      CALL STSetType(st,STSINVERT,IERR)
      CALL STSetShift(st,SHIFT2, IERR)
      CALL STGetKSP(st,ksp,IERR)
!      call KSPCreate(MPI_COMM_WORLD,ksp,ierr)
!      call KSPSetType(ksp,KSPGMRES,ierr)  ! gmres por default
!      call KSPSetFromOptions(ksp,ierr)    ! -ksp_type <name>
!      CALL STSetKSP(st,ksp,IERR)
!      call KSPGetType(ksp,ksptype,ierr)
!      write(*,*) ksptype
      CALL EPSSetFromOptions(eps, IERR)
      call EPSGetProblemType(eps,epsproblemtype,ierr)
      call EPSGetType(eps,epstype,ierr)
      call KSPGetType(ksp,ksptype,ierr)
      call KSPGetPC(ksp,pc,ierr)
      call PCGetType(pc,pctype,ierr)
      call STGetType(st,sttype,ierr)
      call STGetShift(st,shift2,ierr)
      shift=dsqrt(dabs(shift2))
!
      if(irank.eq.0)then
        write(*,'(a)')'************************************************'
        if(epsproblemtype.eq.2)then
          write(*,'(a)')CHAR(9)//'EPS PROBLEM: '//CHAR(9)//'EPS_GHEP'
        else if(epsproblemtype.eq.4)then
          write(*,'(a)')CHAR(9)//'EPS PROBLEM: '//CHAR(9)//'EPS_GNHEP'
        else
          write(*,'(a,i1)') CHAR(9)//'EPS PROBLEM: '//CHAR(9),
     &                      epsproblemtype
        endif
        write(*,'(a)') CHAR(9)//'EPS: '//CHAR(9)//trim(epstype)
        write(*,'(a)') CHAR(9)//'KSP: '//CHAR(9)//trim(ksptype)
        write(*,'(a)') CHAR(9)//'PC: ' //CHAR(9)//trim(pctype)
        write(*,'(a)') CHAR(9)//'ST: ' //CHAR(9)//trim(sttype)
        write(*,'(a,f15.0)') CHAR(9)//'SHIFT: '  //CHAR(9),shift
        write(*,'(a,f15.0)') CHAR(9)//'SHIFT**2: ',shift2
        write(*,'(a)')'************************************************'
      endif
!
      CALL EPSSetBalance(eps,EPS_BALANCE_TWOSIDE,PETSC_DEFAULT_INTEGER,
     &                   PETSC_DEFAULT_REAL,ierr)
      CALL EPSGetBalance(eps,eps_bal_type,its_bal,cutoff,ierr)
      WRITE(EMSG08,*)'-----------------------------------------------'
      CALL PetscPrintf(PETSC_COMM_WORLD,EMSG08//CHAR(10),IERR)
      WRITE(EMSG08,*)' Balancing method: ',eps_bal_type
      CALL PetscPrintf(PETSC_COMM_WORLD,EMSG08//CHAR(10),IERR)
      WRITE(EMSG08,*)' Iterations of the balancing algorithm: ',its_bal
      CALL PetscPrintf(PETSC_COMM_WORLD,EMSG08//CHAR(10),IERR)
      WRITE (EMSG08,*)' Cutoff value: ',cutoff
      CALL PetscPrintf(PETSC_COMM_WORLD,EMSG08//CHAR(10),IERR)
      WRITE (EMSG08,*)'-----------------------------------------------'
      CALL PetscPrintf(PETSC_COMM_WORLD,EMSG08//CHAR(10),IERR)
!
      CALL EPSSolve(eps, IERR)
!
      CALL EPSGetIterationNumber(eps, ITS, IERR)
      CALL EPSGetTolerances(eps, TOL, EPSMAXIT, IERR)
      CALL EPSGetConverged(eps, NCONV, IERR)
      WRITE(ERRMSG,*) ' ITS:',ITS
      CALL PetscPrintf(PETSC_COMM_WORLD, ERRMSG//CHAR(10), IERR)
      WRITE(ERRMSG,*) ' TOL:',TOL
      CALL PetscPrintf(PETSC_COMM_WORLD, ERRMSG//CHAR(10), IERR)
      WRITE(ERRMSG,*) ' EPSMAXIT:',EPSMAXIT
      CALL PetscPrintf(PETSC_COMM_WORLD, ERRMSG//CHAR(10), IERR)
      WRITE(ERRMSG,*) ' NCONV:',NCONV
      CALL PetscPrintf(PETSC_COMM_WORLD, ERRMSG//CHAR(10), IERR)
      TTIME1 = MPI_Wtime()
      WRITE (EMSG04,'(A12,F9.2)') ' TOTAL TIME:', TTIME1 - TTIME0
      CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)
!      
      IF (NCONV .EQ. 0) THEN
         CALL PetscPrintf(PETSC_COMM_WORLD,
     &      ' WARNING: EIGENVALUE NOT CONVERGED!'//CHAR(10), IERR)
         NCONV = 1
      END IF
!      
      open (32,file='autovalores')
      WRITE (32,*)  'FREQUENCY - EIGENVALUE - ERROR'
      WRITE (EMSG08,*)'------------------------------------------------'
      CALL PetscPrintf(PETSC_COMM_WORLD,CHAR(10)//EMSG08//CHAR(10),IERR)
      WRITE (ERRMSG,*)  '** NEV ** FRECUENCIA ** AUTOVALOR ** ERROR **'
      CALL PetscPrintf(PETSC_COMM_WORLD, ERRMSG//CHAR(10), IERR)
      DO ICONT = 0, NCONV - 1
         CALL EPSGetEigenvalue(eps,ICONT,EIG_R,EIG_I,IERR)
       CALL EPSComputeRelativeError(eps,ICONT,EIG_ERROR,IERR)
         WRITE (EMSG12,*) ICONT,dsqrt(EIG_R),EIG_R,
     &                    '+',EIG_I,'i',EIG_ERROR
         CALL PetscPrintf(PETSC_COMM_WORLD, EMSG12//CHAR(10), IERR)
         if(irank.eq.0)then
           WRITE (32,*)  dsqrt(EIG_R),EIG_R,'+',EIG_I,'i',EIG_ERROR
         endif
      ENDDO
      CALL PetscPrintf(PETSC_COMM_WORLD,EMSG08//CHAR(10)//CHAR(10),IERR)
      close(32)
!!!
      DO ICONT = 0, NCONV - 1
         CALL EPSGetEigenpair(eps, ICONT, EIG_R, EIG_I, X, B, IERR)
         CALL EPSComputeRelativeError(eps, ICONT, EIG_ERROR, IERR)
         WRITE (EMSG12,*)  'EIGENVALUE, ERROR:',
     &          EIG_R, '+', EIG_I, 'I, ', EIG_ERROR
         CALL PetscPrintf(PETSC_COMM_WORLD, EMSG12//CHAR(10), IERR)
         WRITE (EMSG08,*)  'FREQUENCY:',dsqrt(EIG_R)
         CALL PetscPrintf(PETSC_COMM_WORLD, EMSG08//CHAR(10), IERR)
         INSERT = 1
         CALL VecScatterBegin(PAR2SEQ, X, IOVEC, INSERT_VALUES, 
     &                                    SCATTER_FORWARD, IERR)
         CALL VecScatterEnd(PAR2SEQ, X, IOVEC, INSERT_VALUES,
     &                                    SCATTER_FORWARD, IERR)
         IF (IRANK .EQ. 0) THEN
           CALL VecGetArray(IOVEC,IOVECSCAN,I_X,IERR)
           CALL GETSOL(IOVECSCAN,I_X,UP,US,UF,PRES,PRESSP,NGLP,
     &                 NGLS,NGLF,NCACOF,nvercos,vec_iperm)
           CALL VecRestoreArray(IOVEC,IOVECSCAN,I_X,IERR)

!*
!******   CALCULO DE LAS SOLUCIONES NODALES EN LA PLACA Y
!*             CREACION DE FICHEROS PARA GRAFICAION.
!*
            IF (INDPLA.EQ.1) THEN   
!               up_max=dabs(up(cnump(nnp(1,1))+1))
!               do i=1,nelp
!                 do j=1,3
!                   if(dabs(up(cnump(nnp(j,i))+1)).gt.up_max)then
!                     up_max=dabs(up(cnump(nnp(j,i))+1))
!                   endif
!                 enddo
!               enddo         
!               up_max=upmax/1000
!               do i=1,nglp
!                 up(i)=up(i)/up_max
!               enddo         
!               do i=1,ngls
!                 us(i)=us(i)/up_max
!               enddo         
!               do i=1,nglf
!                 uf(i)=uf(i)/up_max
!               enddo         
!               do i=1,ncacof
!                 pres(i)=pres(i)/up_max
!               enddo         
!               do i=1,nvercos
!                 pressp(i)=pressp(i)/up_max
!               enddo         
!              CALL vwmpl(NELP,NVERP,NARP,MMP,NNP,CNUMP,NGP,
!     &                     ZP,UP,UPV,BETA_L,ROT,SMP,THIN,ICONT)
              CALL vtkpl(NELP,NVERP,NARP,MMP,NNP,CNUMP,NGP,
     &                     ZP,UP,UPV,BETA_L,ROT,SMP,THIN,ICONT)
              if(ICONT.EQ.0) then
                call vtkpl_case(nconv)
!               call cal_componentesp(nelp,mmp,zp,upv,beta_l,
!     &                              componentes_pz,componentes_prot)
              endif
            ENDIF

!* 
!******   CALCULO DE LAS SOLUCIONES NODALES EN EL SOLIDO Y
!*             CREACION DE FICHEROS PARA GRAFICAION.
!*
            IF (INDSOL.EQ.1) THEN   
              CALL vtkso(NVERS,US,USV,ICONT)
!              CALL vwmso(NVERS,US,USV,ICONT)
              if(ICONT.EQ.0) then
                call vtkso_case(nconv)
              endif
            ENDIF

!*
!******   CALCULO DE LAS SOLUCIONES NODALES EN EL FLUIDO Y
!*             CREACION DE FICHEROS PARA GRAFICAION.
!*
            IF (INDFLU.EQ.1) THEN   
             CALL DEPRE(UF,NELF,NNODF,NVERF,ZF,NGF,MMF,NNF,NSDF,
     &              UFV,pf,IVER,DIV,IOPROM,NDSDF,ISDF,DENSF,VELSQ,VOLF)
               
             if(ICONT.EQ.0) then
               call vtkfl_case(nconv)
!               call desplazamientofluido_vtk_case(nconv)
!              call errestf(nelf,nverf,mmf,ietf,jetf,jeef,zf,nnf,uf,
!     &            ngf,nuf,componentes_f,zp,mmp,nrcof,irefcf,nrcf,ncacof,
!     &            inodfp,numitap,x1p,y1p,x2p,y2p,x3p,y3p,zintp,indcp,
!     &            ielp,dirp,smp,componentes_pz,componentes_prot,ndirf,
!     &            ethau,jlf,uhf,ethaf,EIG_R,velsq,densf,volf,ethaf_aux,
!     &            iver,ethaf_sum)
*              CALL vwmerrfl(nelf,mmf,zf)
!              call solp_noconf(nelf,nglf,mmf,nnf,ngf,zf,volf,densf,velsq
!     &                  ,EIG_R,iflagf,componentes_f,p_noconf,zf_noconf)
!              call solp_conf(nelf,nverf,mmf,nnf,zf,p_noconf,zf_noconf,
!     &                                      volf,vfsum,iver,p_conf)
!              call co_rayleigh(nelf,mmf,zf,p_conf,volf,velsq,rayleigh) !presion
!             call co_rayleigh2(nelf,mmf,zf,ufv,volf,velsq,rayleigh) !desplazamiento
!              open (31,file='output')
!              write(31,310)
! 310          format(1x,'DOF fluido',9x,'Frecuencia',6x,
!     &                  'Cociente de Rayleigh')
!              write(31,311) nglf,dsqrt(EIG_R),dsqrt(rayleigh)
! 311          format(1x,I10,5x,d16.9,5x,d16.9)
             endif
             
*             call desplazamientofluido_vtk(nelf,mmf,zf,nnf,uf,ngf,nuf,
*     &            componentes_f,volf,div,ICONT)
             call vtkfl(nverf,uf,ufv,pf,p_conf,ethaf_aux,ICONT)
            ENDIF

         END IF      !(IRANK.EQ.0)
      ENDDO
!*
!*********************************************************************
!*                  FREEING AND DESTROYING OBJETCS
!*********************************************************************
!*
      CALL VecScatterDestroy(PAR2SEQ, IERR)
      CALL VecDestroy(B, IERR)
      CALL VecDestroy(X, IERR)
      CALL VecDestroy(IOVEC, IERR)
      CALL MatDestroy(KR, IERR) 
      CALL MatDestroy(M, IERR) 
      CALL KSPDestroy(ksp,IERR) 

!*********************************************************************
!     De-Allocamiento de memoria
!*********************************************************************
      DEALLOCATE ( volf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( vfsum, STAT = DEALLOCATEStatus)
      DEALLOCATE ( nuf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( zf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( div, STAT = DEALLOCATEStatus)
      DEALLOCATE ( pf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( ufv, STAT = DEALLOCATEStatus)
      DEALLOCATE ( uf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( iver, STAT = DEALLOCATEStatus)
      DEALLOCATE ( componentes_f, STAT = DEALLOCATEStatus)
      DEALLOCATE ( p_conf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( mmf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( nrcf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( nnf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( nsdf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( ngf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( iflagf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( face2el, STAT = DEALLOCATEStatus)
      DEALLOCATE ( ietf2, STAT = DEALLOCATEStatus)
      DEALLOCATE ( ietf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( jetf, STAT = DEALLOCATEStatus)
      DEALLOCATE ( jeef, STAT = DEALLOCATEStatus)
      DEALLOCATE ( vec_perm, STAT = DEALLOCATEStatus)
      DEALLOCATE ( vec_iperm, STAT = DEALLOCATEStatus)
      IF(DEALLOCATEStatus /= 0) STOP "ERROR in deallocation: "// 
     &  "volf,vfsum,nuf,zf,div,pf,ufv,uf,iver,componentes_f,p_conf,"//
     &  "mmf,nrcf,nnf,nsdf,ngf,iflagf,face2el,ietf2,ietf,"//
     &  "jetf,jeef,vec_perm,vec_iperm"
!*********************************************************************




      TTIME1 = MPI_WTIME()
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' ------------------------------------'//CHAR(10), IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' PAR-IFE3D: NORMAL COMPLETION '//CHAR(10), IERR)
      CALL PetscPrintf(PETSC_COMM_WORLD,
     &     ' ------------------------------------'//CHAR(10), IERR)
      WRITE (EMSG04,'(A12,F9.2)') ' TOTAL TIME:', TTIME1 - TTIME0
      CALL PetscPrintf(PETSC_COMM_WORLD, EMSG04//CHAR(10), IERR)

      CALL SlepcFinalize(IERR)

      STOP
      END
