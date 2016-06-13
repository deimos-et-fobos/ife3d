C
      program ordsur
C
C     Condiciones de contorno
C
C     preprocesador entra los datos de las mallas de vwm
C     y saca mallas para el prog de rt2d (Rodolfo Rodriguez)
C
C     By Ronia-14/3/97
C
      IMPLICIT double precision (a-h,o-z)
      integer elmax
      parameter (nodmax=10000000)
      parameter (maxgroups=100)
      parameter (elmax=1000000)
      
      CHARACTER*80 SUR_OLD,SUR_NEW,SAVESUR,dummy
      dimension mmsur_o(4,elmax),mmsur_n(4,elmax),iaux(2,elmax),
     &          bar_n(3,elmax)
      integer sur_gr_o(maxgroups),sur_gr_n(maxgroups),sur_ngrs,xx,yy
      double precision zf_o(3,nodmax),zf_n(3,nodmax)
      tol=1.d-12

c
c     abre archivo SUR creado por SAVESUR (solo 1 grupo de elementos)
c
      WRITE (6,*) "Archivo sur (SAVESUR): "
      READ  (5,*) SAVESUR
      OPEN  (1,FILE=SAVESUR,STATUS='OLD',ERR=10)

c
c     abre archivo SUR viejo
c
      WRITE (6,*) "Archivo sur (viejo): "
      READ  (5,*) SUR_OLD
      open  (2,FILE=SUR_OLD,STATUS='OLD',ERR=10)

c
c     leo las coordenadas y la conectividad nuevas
c
      call searstr(1,'ELEMENT_GROUPS')
      read(1,*) dummy
      read(1,*) inutil,nelf_n
 
      call searstr(1,'INCIDENCE')
      read(1,*) dummy
      do 1000 i=1,nelf_n
         read(1,*) mmsur_n(1,i),mmsur_n(2,i),mmsur_n(3,i)
         mmsur_n(4,i)=0
 1000 continue
 
      call searstr(1,'COORDINATES')
      read(1,*) nverf_n
      do 1010 i=1,nverf_n
         read(1,*) inutil,zf_n(1,i),zf_n(2,i),zf_n(3,i)
 1010 continue

c
c     leo las coordenadas y la conectividad viejas
c
      call searstr(2,'COORDINATE')
      read(2,*) nverf_o
      do 2030 i=1,nverf_o
         read(2,*) inutil,zf_o(1,i),zf_o(2,i),zf_o(3,i)
 2030 continue

      call searstr(2,'ELEMENT_GROUPS')
      read(2,*) sur_ngrs
      nelf_o=0
      do 2000 i=1,sur_ngrs
            read(2,*) inutil,sur_gr_o(i),dummy
            nelf_o=nelf_o+sur_gr_o(i)
 2000 continue
 
      call searstr(2,'INCIDENCE')
      read(2,*) dummy
      iel=1
      do 2010 i=1,sur_ngrs
         do 2020 j=1,sur_gr_o(i)
            read(2,*) mmsur_o(1,iel),mmsur_o(2,iel),mmsur_o(3,iel)
            mmsur_o(4,iel)=i
            iel=iel+1
 2020    continue
 2010 continue
 


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculo de los baricentros de las caras
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 2300 i0=1,nelf_n
         iaux(1,i0)=0
         iaux(2,i0)=0
         do 2301 j0=1,3
            bar_n(j0,i0)=0
            do 2302 k0=1,3
               bar_n(j0,i0)=bar_n(j0,i0)+zf_n(j0,mmsur_n(k0,i0))
 2302       continue
            bar_n(j0,i0)=bar_n(j0,i0)/3.0
 2301    continue
 2300 continue
! 
!      do 3004 i0=1,nelf_o
!         do 3005 j0=1,3
!            bar_o(j0,i0)=0
!            do 3006 k0=1,3
!               bar_o(j0,i0)=bar_o(j0,i0)+zf_o(j0,mmsur_o(k0,i0))
! 3006       continue
!            bar_o(j0,i0)=bar_o(j0,i0)/3.0
! 3005    continue
! 3004 continue
 
      do 3001 i0=1,nelf_n
         do 3002 j0=1,3
            if ( dabs(zf_n(j0,mmsur_n(1,i0))-
     &                        zf_n(j0,mmsur_n(2,i0))) .lt. tol ) then
              if ( dabs(zf_n(j0,mmsur_n(1,i0))-
     &                        zf_n(j0,mmsur_n(3,i0))) .lt. tol ) then
                iaux(2,i0) = j0
                goto 3066
              endif
            endif
 3002    continue
 3066 continue
 3001 continue
 
      do 3004 i0=1,nelf_o
         do 3005 j0=1,3
            if ( dabs(zf_o(j0,mmsur_o(1,i0))-
     &                        zf_o(j0,mmsur_o(2,i0))) .lt. tol ) then
              if ( dabs(zf_o(j0,mmsur_o(1,i0))-
     &                        zf_o(j0,mmsur_o(3,i0))) .lt. tol ) then
                iaux(1,i0) = j0
                goto 3666
              endif
            endif
 3005    continue
 3666 continue
 3004 continue

      do 4000 i0=1,nelf_n
         do 4010 j0=1,nelf_o
            if (iaux(2,i0).eq.iaux(1,j0)) then
               k0 = iaux(2,i0)
               if (dabs(zf_n(k0,mmsur_n(1,i0))-
     &                        zf_o(k0,mmsur_o(1,j0))) .lt. tol) then
                  xx = 1+mod(k0,3)
                  yy = 1+mod(k0+1,3)
                  det=determinante(zf_o(xx,mmsur_o(1,j0)),
     &                  zf_o(yy,mmsur_o(1,j0)),zf_o(xx,mmsur_o(2,j0)),
     &                  zf_o(yy,mmsur_o(2,j0)),zf_o(xx,mmsur_o(3,j0)),
     &                  zf_o(yy,mmsur_o(3,j0)))
                  aux=determinante(zf_o(xx,mmsur_o(1,j0)),
     &                  zf_o(yy,mmsur_o(1,j0)),zf_o(xx,mmsur_o(2,j0)),
     &                  zf_o(yy,mmsur_o(2,j0)),bar_n(xx,i0),
     &                  bar_n(yy,i0))
                  if ( (aux*det) .gt. 0) then
                  aux=determinante(zf_o(xx,mmsur_o(2,j0)),
     &                  zf_o(yy,mmsur_o(2,j0)),zf_o(xx,mmsur_o(3,j0)),
     &                  zf_o(yy,mmsur_o(3,j0)),bar_n(xx,i0),
     &                  bar_n(yy,i0))
                  if ( (aux*det) .gt. 0) then
                  aux=determinante(zf_o(xx,mmsur_o(3,j0)),
     &                  zf_o(yy,mmsur_o(3,j0)),zf_o(xx,mmsur_o(1,j0)),
     &                  zf_o(yy,mmsur_o(1,j0)),bar_n(xx,i0),
     &                  bar_n(yy,i0))
                  if ( (aux*det) .gt. 0) then
                     mmsur_n(4,i0)=mmsur_o(4,j0)
                     goto 4666
                  endif
                  endif
                  endif
               endif
            endif
 4010    continue
 
 4666    continue
 4000 continue

      do 5000 i0=1,sur_ngrs
         sur_gr_n(i0) = 0
 5000 continue

      do 5010 i0=1,nelf_n
         sur_gr_n(mmsur_n(4,i0)) = sur_gr_n(mmsur_n(4,i0)) + 1
         write(77,*) mmsur_n(4,i0)
 5010 continue

c
c     abre archivo SUR viejo
c
      WRITE (6,*) "Archivo sur (nuevo): "
      READ  (5,*) SUR_NEW
      open  (3,FILE=SUR_NEW,STATUS='UNKNOWN',ERR=10)
      
      write(3,'(a)') '*COORDINATES' 
      write(3,*) nverf_n
      do 6000 i=1,nverf_n
         write(3,*) i,(zf_n(j,i),j=1,3)
 6000 continue

      write(3,'(a)') '*ELEMENT_GROUPS'
      write(3,*) sur_ngrs
      do 6010 i=1,sur_ngrs      
         write(3,123) i,sur_gr_n(i),'Tri3'
 6010 continue
 123  format(1x,i2,1x,i6,1x, a4)
      write(3,'(a)') '*INCIDENCE' 
      write(3,'(a)') '<NONE>'
      nel = 0
      do 6020 i=1,sur_ngrs
         do 6030 j=1,nelf_n
            if(mmsur_n(4,j).eq.i) then
               write(3,*) (mmsur_n(k,j),k=1,3)
               nel = nel + 1
            endif
 6030    continue
 6020 continue
      
      if(nel.ne.nelf_n) then
         write (6,*) "Error. nel = ",nel,"; nelf_n = ",nelf_n
         STOP
      endif

      write(3,'(a)') '*END'
      close(3)
      close(1)
      close(2)
      return

 10   write (6,*) "Error al abrir el archivo de preprocesador"
      STOP
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                        FUNCION DETERMINANTE	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function determinante(x1,y1,x2,y2,x3,y3)
      implicit double precision (a-h,o-z) 
      determinante=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      end
