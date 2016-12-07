!************************************************************************
!*
!*                   LECTURA DE LAS MALLAS
!*
!************************************************************************

      subroutine lecmalla (mallap,mallas,mallaf,indpla,indsol,indflu,
     &                dirp,smp,thin,nelp,narp,nverp,mmp,nnp,nrap,nrep,
     &                  nrvp,zp,nsdp,nels,nvers,mms,nrvs,nrcs,zs,nsds,
     &                  nelf,nnodf,nverf,nnf,mmf,nrcf,zf,nsdf)

      implicit double precision (a-h,o-z)

      dimension mmp(3,*),nnp(6,*),nrep(*),nrap(3,*),nrvp(3,*),nsdp(*),
     &          zp(2,*),thin(*),smp(*),mms(4,*),nrvs(4,*),nrcs(4,*),
     &          nsds(*),zs(3,*),mmf(4,*),nrcf(4,*),nnf(4,*),zf(3,*),
     &          nsdf(*)
      integer dirp(*)
      
      character*80 mallap,mallas,mallaf,none

!*
!******  Lectura de la malla de la placa.
!*
      if (indpla.eq.0) then
         nelp=0
         nverp=0
         narp=0
      else
         open(10,file=mallap)
         write(6,*) "---> ",'Reading PLATE mesh'
         call searstr(10, 'COORDINATES')
         read(10,*,err=3020,end=3020) nverp
         do i=1,nverp
            read(10,*,err=3020,end=3020)   j,zp(1,i),zp(2,i)
         enddo
         nelp=0
         call searstr(10, 'ELEMENT_GROUPS')
         read(10,*)ngroupp
         do 1000 i=1,ngroupp
            read(10,*) aux, nel_gr, none
            nelp=nelp+nel_gr
 1000    continue
         call searstr(10, 'INCIDENCE')
         read(10,'(a)') none
         read(10,*,err=3020,end=3020) (( mmp(i,k),i=1,3),k=1,nelp)
         call searstr(10, 'NNP')
         read(10,*) nverp
         read(10,*) narp
         call searstr(10, 'NODES_INCIDENCE')
         read(10,*,err=3020,end=3020) (( nnp(i,k),i=1,6),k=1,nelp)
         call searstr(10, 'NRE')
         do 1001 k=1,nelp
           read(10,*,err=3020,end=3020) nrep(k),nsdp(k),dirp(k),
     &                                  smp(k),thin(k)
 1001    continue
!         call searstr(10, 'NRE')
!         read(10,*,err=3020,end=3020) (nrep(k),k=1,nelp)
!         call searstr(10, 'NSD')
!         read(10,*,err=3020,end=3020) ( nsdp(k),k=1,nelp)
!         call searstr(10, 'DIRP')
!         read(10,*,err=3020,end=3020) (dirp(k),k=1,nelp)
!         call searstr(10, 'SMP')
!         read(10,*,err=3020,end=3020) (smp(k),k=1,nelp)
!         call searstr(10, 'THIN')
!         read(10,*,err=3020,end=3020) (thin(k),k=1,nelp)
         call searstr(10, 'NRA')
         read(10,*,err=3020,end=3020) ((nrap(i,k),i=1,3),k=1,nelp)
         call searstr(10, 'NRV')
         read(10,*,err=3020,end=3020) ((nrvp(i,k),i=1,3),k=1,nelp)
         close(10)
      endif   

!*
!******  Lectura de la malla del sólido.
!*
      if (indsol.eq.0) then
         nels=0
         nvers=0
      else
         open(11,file=mallas)
         write(6,*) "---> ",'Reading SOLID mesh'
         call searstr(11, 'COORDINATES')
         read(11,*,err=3020,end=3020) nvers
         do i=1,nvers
            read(11,*,err=3020,end=3020)   j,zs(1,i),zs(2,i),zs(3,i)
         enddo
         nels=0
         call searstr(11, 'ELEMENT_GROUPS')
         read(11,*)ngroups
         do 1010 i=1,ngroups
            read(11,*) aux, nel_gr, none
            nels=nels+nel_gr
 1010    continue
         call searstr(11, 'INCIDENCE')
         read(11,'(a)') none
         read(11,*,err=3020,end=3020) (( mms(i,k),i=1,4),k=1,nels)
         call searstr(11, 'NRV')
         read(11,*,err=3020,end=3020) ((nrvs(i,k),i=1,4),k=1,nels)
         call searstr(11, 'NRC')
         read(11,*,err=3020,end=3020) ((nrcs(i,k),i=1,4),k=1,nels)
         call searstr(11, 'NSD')
         read(11,*,err=3020,end=3020) ( nsds(k),k=1,nels)
         close(11)
      endif 

!*
!******  Lectura de la malla del fluido.
!*
      if (indflu.eq.0) then
         nelf=0
         nnodf=0
         nverf=0
      else
         open(12,file=mallaf)
         write(6,*) "---> ",'Reading FLUID mesh'
         call searstr(12, 'COORDINATES')
         read(12,*,err=3020,end=3020) nverf
         do 1050 i=1,nverf
            read(12,*,err=3020,end=3020)   j,zf(1,i),zf(2,i),zf(3,i)
 1050    continue
         nelf=0
         call searstr(12, 'ELEMENT_GROUPS')
         read(12,*)ngroupf
         do 1100 i=1,ngroupf
            read(12,*) aux, nel_gr, none
            nelf=nelf+nel_gr
 1100    continue
         call searstr(12, 'INCIDENCE')
         read(12,'(a)') none
         read(12,*,err=3020,end=3020) (( mmf(i,k),i=1,4),k=1,nelf)
         call searstr(12, 'NFACES')
         read(12,*,err=3020,end=3020) nnodf
         call searstr(12, 'FACE_INCIDENCE')
         read(12,*,err=3020,end=3020)(( nnf(j,k),j=1,4),k=1,nelf)
         call searstr(12, 'NRC')
         read(12,*,err=3020,end=3020)((nrcf(j,k),j=1,4),k=1,nelf)
         call searstr(12, 'NSD')
         read(12,*,err=3020,end=3020)(nsdf(k),k=1,nelf)
         close(12)
      endif

!*
!******  Fin de lectura de mallas
!*
      return

 3020 WRITE (6, 2020)
 2020 FORMAT (//,' Error in data file')
      STOP ' '

      end
