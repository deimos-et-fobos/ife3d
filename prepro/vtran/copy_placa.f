
      subroutine copy_placa(mm,mmsur,z,ncopy,dx,dy,dz)

      IMPLICIT double precision (a-h,o-z)
      parameter (maxgroups=1000)
      
      CHARACTER*80 FILEBC,FILEVWM,FILESUR,VWMNEW,SURNEW,BCNEW
      CHARACTER*180 dummy1,dummy2
      dimension mm(4,*),mmsur(3,*),z(3,*),pos(2,maxgroups)
      integer el_gr(maxgroups),sur_gr(maxgroups),el_ngrs,sur_ngrs
      integer  ncopy,NR(maxgroups),elemref(3,maxgroups)

!     abre archivo VWM
      WRITE (6,*) "Archivo vwm: "
      READ  (5,*) FILEVWM
      OPEN  (11,FILE=FILEVWM,STATUS='OLD',ERR=10)

!     abre archivo SUR
      WRITE (6,*) "Archivo sur: "
      READ  (5,*) FILESUR
      open  (12,FILE=FILESUR,STATUS='OLD',ERR=10)

!     abre archivo BC 
      WRITE (6,*) "¿Modificar el archivo bc? (si=1;no=0)"
      READ  (5,*) nbc
      if(nbc.eq.1)then
        WRITE (6,*) "Archivo bc: "
        READ  (5,*) FILEBC
        OPEN  (13,FILE=FILEBC,STATUS='OLD',ERR=10)
      endif

      call searstr(11,'COORDINATES')
      read(11,*) nver
      do 100 i=1,nver
        read(11,*) inutil,z(1,i),z(2,i)
 100  continue

      call searstr(11,'ELEMENT_GROUPS')
      read(11,*) el_ngrs
      nel=0
      do 200 i=1,el_ngrs
        read(11,*) inutil, el_gr(i)
        nel=nel+el_gr(i)
 200  continue

      call searstr(11,'INCIDENCE')
      read(11,*)
      do 300 i=1,nel
        read(11,*) mm(1,i),mm(2,i),mm(3,i)
 300  continue

      call searstr(12,'ELEMENT_GROUPS')
      read(12,*)sur_ngrs
      nelsur=0
      do 400 i=1,sur_ngrs
        read(12,*) inutil, sur_gr(i)
        nelsur=nelsur+sur_gr(i)
 400  continue  

      call searstr(12,'INCIDENCE')
      read(12,*)
      do 500 i=1,nelsur
        read(12,*) mmsur(1,i),mmsur(2,i)
 500  continue

      if(nbc.eq.1)then
        read(13,'(a180)',end=20) dummy1
        do 600 i=1,sur_ngrs
          read(13,*,end=20) NG, NR(i) 
 600    continue   
        read(13,*,end=20)
        read(13,'(a180)',end=20) dummy2
        do 610 i=1,el_ngrs
          read(13,*,end=20) elemref(1,i),elemref(2,i),elemref(3,i),
     &                      pos(1,i),pos(2,i)
 610    continue
      endif

!     Escritura de archivos nuevos
      VWMNEW = trim(FILEVWM) // '.new'
      SURNEW = trim(FILESUR) // '.new'
      open(14,file=VWMNEW)
      open(15,file=SURNEW)
      if(nbc.eq.1)then
        BCNEW = trim(FILEBC) // '.new'
        open(16,file=BCNEW)
      endif

      write(14,'(a)') '*COORDINATES'
      write(15,'(a)') '*COORDINATES'
      write(14,*) nver*ncopy
      write(15,*) nver*ncopy
      do 1000 k=0,ncopy-1 
        do 1100 i=1,nver
          write(14,*) k*nver+i,z(1,i)+k*dx,z(2,i)+k*dy
          write(15,*) k*nver+i,z(1,i)+k*dx,z(2,i)+k*dy
 1100   continue
 1000 continue
 
      write(14,'(a)') '*ELEMENT_GROUPS'
      write(14,*) el_ngrs*ncopy
      do 2000 k=0,ncopy-1
        do 2100 i=1,el_ngrs
          write(14,123) k*el_ngrs+i,el_gr(i),'Tri3  '
 2100   continue
 2000 continue
 123  format(1x,i3,1x,i9,1x,a6)

      write(14,'(a)') '*INCIDENCES' 
      write(14,'(a)') '<NONE>'
      do 3000 k=0,ncopy-1
        do 3100 i=1,nel
          write(14,*) mm(1,i)+k*nver,mm(2,i)+k*nver,
     &                mm(3,i)+k*nver
 3100   continue
 3000 continue

      write(15,'(a)') '*ELEMENT_GROUPS'
      write(15,*) sur_ngrs*ncopy
      do 4000 k=0,ncopy-1
        do 4100 i=1,sur_ngrs
          write(15,123) k*sur_ngrs+i,sur_gr(i),'Seg2  '
 4100   continue
 4000 continue

      write(15,'(a)') '*INCIDENCES' 
      write(15,'(a)') '<NONE>'
      do 5000 k=0,ncopy-1
        do 5100 i=1,nelsur
          write(15,*) mmsur(1,i)+k*nver,mmsur(2,i)+k*nver
 5100   continue
 5000 continue

      write(14,'(a)') '*END'
      write(15,'(a)') '*END'

      if(nbc.eq.1)then
        write(16,'(a180)') dummy1
        do 6000 k=0,ncopy-1
          do 6100 i=1,sur_ngrs
            write(16,*) i+k*sur_ngrs, NR(i) 
 6100     continue
 6000   continue
        write(16,'(a180)') ' '
        write(16,'(a)') dummy2
        do 6200 k=0,ncopy-1
          do 6300 i=1,el_ngrs
            write(16,*) i+k*el_ngrs,elemref(2,i),elemref(3,i),
     &                  pos(1,i)+k*dz,pos(2,i)
 6300     continue
 6200   continue
        close(13)
        close(16)
      endif

      close(11)
      close(12)
      close(14)
      close(15)

      return

 10   write (6,*) "Error al abrir el archivo de preprocesador"
      STOP
 20   write(*,*) "Error al leer los limites de las referencias"
      STOP
      end







