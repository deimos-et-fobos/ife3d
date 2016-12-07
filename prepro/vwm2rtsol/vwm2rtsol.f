C
      program vwm2rtsol
C
C     Condiciones de contorno
C
C     preprocesador
C
C
C     By Ronia-11/3/97
C
      IMPLICIT double precision (A-H,O-Z)

      parameter (maxgroups=1000)
      parameter (maxnod=10000000)
      parameter (maxelem=20000000)
      
      dimension mms(4,maxelem),mmsur(4,maxelem),NR(maxgroups),
     &       zs(3,maxnod),nrvs(4,maxelem),nrcs(4,maxelem),je(4*maxelem),
     &       nref(maxnod),iet(maxnod),jet(4*maxelem),jee(4*maxelem)
      CHARACTER*80 FILESUR,FILEVWM,FILEBO,FILEOUT
      integer el_gr(maxgroups),sur_gr(maxgroups),el_ngrs,sur_ngrs,
     &        xx(8,8)
      data xx/0,1,2,3,4,5,6,7,       1,1,1,1,1,1,1,1,
     &        2,1,2,5,7,5,1,7,       3,1,5,3,6,5,6,1,
     &        4,1,7,6,4,1,6,7,       5,1,5,5,1,5,1,1,
     &        6,1,1,6,6,1,6,1,       7,1,7,1,7,1,1,7/

      write(*,*)"Starting vwm2rtsol ..."
      write(*,*)" Valid BC:"
      write(*,*)"  0: Libre"
      write(*,*)"  1: Ux = Uy = Uz (Empotramiento)"
      write(*,*)"  2: Ux = 0 (Sin desplazamiento en direccion x)"
      write(*,*)"  3: Uy = 0 (Sin desplazamiento en direccion y)"
      write(*,*)"  4: Uz = 0 (Sin desplazamiento en direccion z)"
      write(*,*)"  5: Ux = Uy = 0 " 
      write(*,*)"  6: Uy = Uz = 0 "
      write(*,*)"  7: Uz = Ux = 0 "
      write(*,*)" >7: Contacto"
      write(*,*)

c
c     abre archivo VWM
c
!      WRITE (6,*) "Archivo vwm: "
!      READ  (5,*) FILEVWM
      FILEVWM='solido.vwm'
      OPEN  (1,FILE=FILEVWM,STATUS='OLD',ERR=210)
      write(*,*)"Abriendo solido.vwm ..."

c
c     abre archivo SUR
c
!      WRITE (6,*) "Archivo sur: "
!      READ  (5,*) FILESUR
      FILESUR='solido.sur'
      OPEN  (2,FILE=FILESUR,STATUS='OLD',ERR=220)
      write(*,*)"Abriendo solido.sur ..."
            
c
c     abre archivo de datos
c
!      WRITE (6,*) "Archivo de datos del preprocesador: "
!      READ  (5,*) FILEBO
      FILEBO='bc-so.dat'
      OPEN  (4,FILE=FILEBO,STATUS='OLD',ERR=230)
      write(*,*)"Abriendo bc-so.dat ..."
      
c
c     Scanea el nombre del archivo de salida
c
!      WRITE (6,*) "Archivo salida: "
!      READ  (5,*) FILEOUT
      FILEOUT='solido.rt'
      open(12,file=FILEOUT,ERR=240)
      write(*,*)"Abriendo solido.rt ..."

c
c     abro el archivo con la malla del solido
c
      call searstr(1,'COORDINATES')
      read(1,*) nvers
      do 6020 i=1,nvers
        read(1,*) inutil ,zs(1,i),zs(2,i),zs(3,i)
 6020 continue

      call searstr(1,'ELEMENT_GROUPS')
      read(1,*)el_ngrs
      nels=0
      do 6000 i=1,el_ngrs
        read(1,*) inutil ,el_gr(i)
        nels=nels+el_gr(i)
 6000 continue

c
c     Leo la matriz de conectividad debe ser <NONE>
c
      call searstr(1,'INCIDENCE')
      read(1,*)
      do 6010 i=1,nels
        read(1,*) mms(1,i),mms(2,i),mms(3,i),mms(4,i)
 6010 continue

c
c     abro el archivo con la malla de superficie del solido
c 
      call searstr(2,'ELEMENT_GROUPS')
      read(2,*)sur_ngrs
      nelsur=0
      do 1100 i0=1,sur_ngrs
        read(2,*) inutil, sur_gr(i0)
        nelsur=nelsur+sur_gr(i0)
 1100 continue  

      call searstr(4,'BC_ELEMENT_GROUPS')
      do 1200 i=1,sur_ngrs
        read(4,*,end=200) NG,NR(i)
 1200 continue
    
c
c     Leo la matriz de conectividad de la superficie
c
      call searstr(2,'INCIDENCE')
      read(2,*) 
      iel=1
      do 1300 i=1,sur_ngrs
        do 1400 j=1,sur_gr(i)
          read(2,*) mmsur(1,iel),mmsur(2,iel),mmsur(3,iel)
          mmsur(4,iel)=NR(i)
          iel=iel+1
 1400   continue
 1300 continue
 
      close(1)
      close(2)
      close(4)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c		Inicializa las matrices 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!     pone a cero la matriz de vertices, aristas y vector de ref de elementos
    
      do 6030 i=1,nels
        do 6040 j=1,4
          nrvs(j,i)=0
          nrcs(j,i)=0
          je(4*(i-1)+j)=mms(j,i)
 6040   continue
 6030 continue

!	pone a cero los vectores de referencias de bordes     
      do 6050 i=1, nvers
        nref(i)=0
 6050 continue  
     
      call trasim(4,je,iet,jet,nvers,nels)
      CALL conel4(je,iet,jet,jee,nels)
      
c
c	Calculo de la matriz  nrcs
c
      do 3010 i=1,nelsur
        ielem=0
        N1=mmsur(1,i)
        N2=mmsur(2,i)
        N3=mmsur(3,i)
        do 3020 j=iet(N1),iet(N1+1)-1
          iev = jet(j)
          in = 4*iev-3           
          do 3030 k=in,in+3
            if(je(k).eq.N2)then
              do 3040 l=in,in+3
                if (je(l).eq.N3) then
                  ielem=iev
                endif
 3040         continue
            endif
 3030     continue
 3020   continue
        if (ielem.ne.0) then
          do 3100 j=1,4
            if( mms(j,ielem).ne.N1 .and. mms(j,ielem).ne.N2
     &              .and. mms(j,ielem).ne.N3 )then
              nrcs(j,ielem)=mmsur(4,i)
              goto 3600
            endif
 3100     continue
        endif
 3600   continue
        if (mmsur(4,i).le.7) then
          do 3050 j=1,3
            nref(mmsur(j,i))=xx(nref(mmsur(j,i))+1,mmsur(4,i)+1)
 3050     continue
        endif
 3010 continue

c
c      Recorrer elementos, setear matriz referencias mediante nref
c
      do 5000 i=1,nels
         do 5010 j=1,4
            nrvs(j,i)=nref(mms(j,i))  
 5010    continue
 5000 continue 

c
c     Archivo de salida de datos
c

      write(12,'(a)') '*COORDINATES' ! Coordenadas
      write(12,*) nvers
      do 4000 i=1,nvers
         write(12,*) i,zs(1,i),zs(2,i),zs(3,i)
 4000 continue
      write(12,*)

      write(12,'(a)') '*ELEMENT_GROUPS'
      write(12,*) el_ngrs
      do 4010 i=1,el_ngrs
            write(12,123) i,el_gr(i),'Tetra4'
 4010 continue
      write(12,*)
 123  format(1x,i1,1x,i10,1x,a6)
      write(12,'(a)') '*INCIDENCE' ! Coordenadas
      write(12,'(a)') '<NONE>'
      do 4020 i=1,nels
         write(12,*) mms(1,i),mms(2,i),mms(3,i),mms(4,i)
 4020 continue
      write(12,*)

      write(12,'(a)') '*NRV'    
      do 4030 i=1,nels
         write(12,*) nrvs(1,i),nrvs(2,i),nrvs(3,i),nrvs(4,i)
 4030 continue
      write(12,*)

      write(12,'(a)') '*NRC'    
      do 4040 i=1,nels
         write(12,*) nrcs(1,i),nrcs(2,i),nrcs(3,i),nrcs(4,i)
 4040 continue
      write(12,*)
         
      write(12,'(a)') '*NSD'   !Numero de Subdominios
      write(12,*)  (' 1 ',i=1,nels)
      write(12,*)

      write(12,'(a)') '*END' ! Final
      close(12)

      return

 210  write (6,*) "Error al abrir el archivo: 'solido.vwm'"
      STOP
 220  write (6,*) "Error al abrir el archivo: 'solido.sur'"
      STOP
 230  write (6,*) "Error al abrir el archivo: 'bc-so.dat'"
      STOP
 240  write (6,*) "Error al abrir el archivo: 'solido.rt'"
      STOP
 200  write(*,*) "Error al leer los limites de las referencias"
      STOP
      end
