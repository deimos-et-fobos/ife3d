C
      program vwm2rtpla
C
C     Condiciones de contorno
C
C     preprocesador
C
C
C     By Ronia-11/3/97
C
      IMPLICIT double precision (A-H,O-Z)
      
      parameter (maxlist=1300)
      parameter (maxgroups=1000)
      parameter (maxnod=10000000)
      parameter (maxelem=10000000)
      
      dimension mmp(3,maxelem),mmsur(3,maxelem),NR(maxgroups),nn(4),
     &     nnp(6,maxelem),zp(2,maxnod),nrvp(3,maxelem),nrap(3,maxelem),
     &     nrep(maxelem),nref(maxnod),smp(maxgroups),thin(maxgroups),
     &     smpp(maxelem),thinp(maxelem),nsd(maxgroups),nsdp(maxelem),
     &     je(3*maxnod),iet(maxnod),jet(3*maxnod),jee(3*maxelem)
      CHARACTER*80 FILESUR,FILEVWM,FILEBO,FILEOUT,dummy 
      integer el_gr(maxgroups),sur_gr(maxgroups),el_ngrp,sur_ngrp,
     &        dirp(maxgroups),dirps(maxelem),xx(9,9)
      data xx/0,1,2,3,4,5,6,7,8,       1,1,1,1,1,1,1,1,1,
     &        2,1,2,9,7,11,11,7,2,     3,1,9,3,10,6,6,10,3,
     &        4,1,7,10,4,12,1,7,7,     5,1,11,6,12,5,6,1,6,
     &        6,1,11,6,1,6,6,1,6,      7,1,7,10,7,1,1,7,7,
     &        8,1,2,3,7,6,6,7,8/

      write(*,*)"Starting vwm2rtpla ..."
      write(*,*)" Valid BC:"
      write(*,*)"  0: Libre"
      write(*,*)"  1: Ux, Uy, Uz, Rx, Ry = 0 (Empotramiento)"
      write(*,*)"  2: Uz, Rx = 0 (Simplemente apoyado en direccion x)"
      write(*,*)"  3: Uz, Ry = 0 (Simplemente apoyado en direccion y)"
      write(*,*)"  4: Ux, Rx = 0 (Simetría respecto al eje y)"
      write(*,*)"  5: Uy, Ry = 0 (Simetría respecto al eje x)"
      write(*,*)"  6: Uz, Uy, Ry = 0 (Antisimetría respecto al eje y)"
      write(*,*)"  7: Uz, Ux, Rx = 0 (Antisimetría respecto al eje x)"
      write(*,*)"  8: Uz = 0 (Simplemente apoyado sin restingir Rx,Ry)"
      write(*,*)">12: Contacto"
      write(*,*)

c
c	abre archivo VWM
c
!      WRITE (6,*) "Archivo vwm: "
!      READ  (5,*) FILEVWM
      FILEVWM='placa.vwm'
      OPEN  (1,FILE=FILEVWM,STATUS='OLD',ERR=210)
      write(*,*)"Abriendo placa.vwm ..."

c
c	abre archivo SUR
c
!      WRITE (6,*) "Archivo sur: "
!      READ  (5,*) FILESUR
      FILESUR='placa.sur'
      open  (2,FILE=FILESUR,STATUS='OLD',ERR=220) 
      write(*,*)"Abriendo placa.sur ..."
      

c
c	abre archivo de datos
c
!      WRITE (6,*) "Archivo de datos del preprocesador: "
!      READ  (5,*) FILEBO
      FILEBO='bc-pl.dat'
      OPEN  (4,FILE=FILEBO,STATUS='OLD',ERR=230)
      write(*,*)"Abriendo bc-pl.dat ..."
      
c
c	Scanea el nombre del archivo de salida
c
!      WRITE (6,*) "Archivo salida: "
!      READ  (5,*) FILEOUT
      FILEOUT='placa.rt'
      open(12,file=FILEOUT,ERR=240)
      write(*,*)"Abriendo placa.rt ..."

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Lee el archivo VWM
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c
c     abro el archivo con las malla del solido
c

      call searstr(1,'ELEMENT_GROUPS')
      read(1,*)el_ngrp
      nelp=0
      do 6000 i=1,el_ngrp
            read(1,*) inutil ,el_gr(i)
            nelp=nelp+el_gr(i)
 6000 continue
 
c
c     Leo la matriz de conectividad debe ser <NONE>
c
      call searstr(1,'INCIDENCES')
      read(1,*)
      do 6010 i=1,nelp
         read(1,*) mmp(1,i),mmp(2,i),mmp(3,i)
 6010 continue
c
c     lee las coordenadas
c
  
      call searstr(1,'COORDINATES')
      read(1,*) nverp
      do 6020 i=1,nverp
         read(1,*) inutil ,zp(1,i),zp(2,i)
 6020 continue

c
c     abro el archivo con las malla del borde del solido
c 
      call searstr(2,'ELEMENT_GROUPS')
      read(2,*)sur_ngrp
      nelsur=0
      do 1100 i0=1,sur_ngrp
            read(2,*) inutil, sur_gr(i0)
            nelsur=nelsur+sur_gr(i0)
 1100 continue  

      call searstr(4,'BC_ELEMENT_GROUPS')
      do 1200 i=1,sur_ngrp
        read(4,*,end=200) NG,NR(i)
 1200 continue   
c
c     Leo la matriz de conectividad del borde
c
      call searstr(2,'INCIDENCES')
      read(2,*) dummy
      iel=1
      do 1300 i=1,sur_ngrp
            do 1400 j=1,sur_gr(i)
                  read(2,*) mmsur(1,iel),mmsur(2,iel)
                  mmsur(3,iel)=NR(i)
                  iel=iel+1
 1400       continue
 1300 continue
      
c
c     Calculo nnp
c
      do 2020 i=1,nelp
         do 2010 j=1,3
            nnp(j,i)=mmp(j,i)
            nnp(j+3,i)=0
            je(3*(i-1)+j)=mmp(j,i)
 2010    continue
 2020 continue
 
      CALL trasim(3,je,iet,jet,nverp,nelp)
      CALL conel3(je,iet,jet,jee,nelp)
      
      icont=nverp
      do 2030 i=1,nelp
        do 2040 j=1,3
          if(jee(3*(i-1)+j).eq.0)then
            icont=icont+1
            nnp(j+3,i)=icont
          else
            if(jee(3*(i-1)+j).gt.i)then
              icont=icont+1
              ivec=jee(3*(i-1)+j)
              nnp(j+3,i)=icont
              do 2050 k=1,3
                if(jee(3*(ivec-1)+k).eq.i) nnp(k+3,ivec)=icont
 2050         continue
            endif
          endif
 2040   continue
 2030 continue
      nnodp=icont
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Inicializa las matrices 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!     pone a cero la matriz de vertices, aristas y vector de ref de elementos
      do 6030 i=1, nelp
         do 6040 j=1,3
            nrvp(j,i)=0
            nrap(j,i)=0
 6040    continue
         nrep(i)=0
 6030 continue
      narp=nnodp-nverp
      
!	pone a cero los vectores de referencias de bordes     
      do 6050 i=1, nverp
         nref(i)=0
 6050 continue  
 
c
c	Calculo de la matriz  nrap
c

      do 1010 i=1,nelp
        nn(1)=mmp(1,i)
        nn(2)=mmp(2,i)
        nn(3)=mmp(3,i)
        nn(4)=nn(1)            
        do 1020 j=1,3
          do 1030 k=1,nelsur
            n1=mmsur(1,k)
            n2=mmsur(2,k)
            if (nn(j).eq.n1) then
              if (nn(j+1).eq.n2) nrap(j,i)=mmsur(3,k)
            end if

            if((mmsur(3,k).ge.0).and.(mmsur(3,k).le.12))then
              nref(n1)=xx(nref(n1)+1,mmsur(3,k)+1)
              nref(n2)=xx(nref(n2)+1,mmsur(3,k)+1)
            endif
!            if (mmsur(3,k).eq.1) then
!              nref(n1)=1
!              nref(n2)=1
!            else
!              if (mmsur(3,k).eq.2) THEN
!                if (nref(n1).eq.0) nref(n1)=2
!                if (nref(n2).eq.0) nref(n2)=2
!                if (nref(n1).eq.3) nref(n1)=1
!                if (nref(n2).eq.3) nref(n2)=1
!                if (nref(n1).eq.4) nref(n1)=2
!                if (nref(n2).eq.4) nref(n2)=2
!                if (nref(n1).eq.5) nref(n1)=1
!                if (nref(n2).eq.5) nref(n2)=1
!              else
!                if (mmsur(3,k).eq.3) THEN
!                  if (nref(n1).eq.0) nref(n1)=3
!                  if (nref(n2).eq.0) nref(n2)=3
!                  if (nref(n1).eq.2) nref(n1)=1
!                  if (nref(n2).eq.2) nref(n2)=1
!                  if (nref(n1).eq.4) nref(n1)=1
!                  if (nref(n2).eq.4) nref(n2)=1
!                  if (nref(n1).eq.5) nref(n1)=3
!                  if (nref(n2).eq.5) nref(n2)=3
!                else
!                  if (nref(n1).eq.0) nref(n1)=mmsur(3,k)
!                  if (nref(n2).eq.0) nref(n2)=mmsur(3,k)
!                endif
!              endif
!            end if  
 1030     continue
 1020   continue
 1010 continue

c
c      Recorrer elementos, setear matriz referencias mediante nref
c
      do 5000 i=1,nelp
         do 5010 j=1,3
            nrvp(j,i)=nref(mmp(j,i))  
 5010    continue
 5000 continue 
 
      iel=1
      call searstr(4,'ELEMENT_GROUPS')
      do 3010 i=1,el_ngrp
        read(4,*,end=200) NG,NR(i),nsd(i),dirp(i),smp(i),thin(i)
        do 3020 i0=1,el_gr(i)
          if (NR(i).NE.0) nrep(iel)=NR(i)
          nsdp(iel)=nsd(i)
          dirps(iel)=dirp(i)
          smpp(iel)=smp(i)
          thinp(iel)=thin(i)
          iel=iel+1
 3020   continue
 3010 continue
 

c
c     Archivo de salida de datos
c

      write(12,'(a)') '*COORDINATES' ! Coordenadas
      write(12,*) nverp
      do 4000 i=1,nverp
         write(12,*) i,zp(1,i),zp(2,i)
 4000 continue
      write(12,*)

      write(12,'(a)') '*ELEMENT_GROUPS'
      write(12,*) el_ngrp
      do 4010 i=1,el_ngrp
            write(12,123) i,el_gr(i),'TRI3'
 4010 continue
      write(12,*)
 123  format(1x,i3,1x,i6,1x, a4)
      write(12,'(a)') '*INCIDENCE' ! Coordenadas
      write(12,'(a)') '<NONE>'
      do 4020 i=1,nelp
         write(12,*) mmp(1,i),mmp(2,i),mmp(3,i)
 4020 continue
      write(12,*)

      write(12,'(a)') '*NNP'
      write(12,*) nverp
      write(12,*) narp
      write(12,*)
 
      write(12,'(a)') '*NODES_INCIDENCE' ! Incidencias de nodos
      do 4030 i=1,nelp
         write(12,*)  (nnp(j,i),j=1,6) 
 4030 continue
      write(12,*)
   
      write(12,'(a)') ' NRE | NSD | DIRP(1=x,2=y,3=z) | SMP | THICKNESS'
      write(12,'(a)') '*NRE'  
      do 4031 i=1,nelp
        write(12,*)  nrep(i),nsdp(i),dirps(i),smpp(i),thinp(i) 
 4031 continue
      write(12,*)

!      write(12,'(a)') '*NRE'  
!      write(12,*)  (nrep(i),i=1,nelp) 
!      write(12,*)
!      
!      write(12,'(a)') '*DIRP'  
!      write(12,*)  (dirps(i),i=1,nelp) 
!      write(12,*)
!      
!      write(12,'(a)') '*SMP'  
!      write(12,*)  (smpp(i),i=1,nelp) 
!      write(12,*)
!      
!      write(12,'(a)') '*THIN'  
!      write(12,*)  (thinp(i),i=1,nelp)  
!      write(12,*)
      
      write(12,'(a)') '*NRV'    
      do 4040 i=1,nelp
         write(12,*) nrvp(1,i),nrvp(2,i),nrvp(3,i)
 4040 continue
      write(12,*)

      write(12,'(a)') '*NRA'    
      do 4050 i=1,nelp
        write(12,*) nrap(1,i),nrap(2,i),nrap(3,i)
 4050 continue
      write(12,*)
         
!      write(12,'(a)') '*NSD'   !Numero de Subdominios
!      write(12,*)  (' 1 ',i=1,nelp)
!      write(12,*)

      write(12,'(a)') '*END' ! Final


c
c     Cierra los archivos de datos
c

      close(1)
      close(2)
      close(3)
      close(4)
      close(12)

      return

 210  write (6,*) "Error al abrir el archivo: 'placa.vwm'"
      STOP
 220  write (6,*) "Error al abrir el archivo: 'placa.sur'"
      STOP
 230  write (6,*) "Error al abrir el archivo: 'bc-pl.dat'"
      STOP
 240  write (6,*) "Error al abrir el archivo: 'placa.rt'"
      STOP
 200  write(*,*) "Error al leer los limites de las referencias"
      STOP
      end
