C
      program vwm2rtflu_mpi
C
C     Condiciones de contorno
C
C     preprocesador entra los datos de las mallas de vwm
C     y saca mallas para el prog de rt2d (Rodolfo Rodriguez)
C
C     By Ronia-14/3/97
C
      IMPLICIT double precision (a-h,o-z)
      include 'mpif.h'
      parameter (maxgroups=1000)
      parameter (maxnod=4000000)
      parameter (maxelem=17000000)
      parameter (maxface=4000000)
      
      CHARACTER*80 FILEBO,FILEVWM,FILESUR,FILEOUT,dummy
      dimension mmf(4,maxelem),nnf(4,maxelem),mmsur(4,maxelem),
     &          zf(3,maxnod),nrcf(4,maxelem),NR(maxgroups),
     &          je(4*maxelem),iet(maxnod),jet(4*maxelem),jee(4*maxelem),
     &          zbar(3,maxface)
      integer cn(4,3),cont,el_gr(maxgroups),sur_gr(maxgroups),
     &            el_ngrs,sur_ngrs,vecino,flag(maxface)
      integer rank, size, mpierr, tag, elemperproc 
      REAL*4 TI,TA(2)
      tol=1.d-6

      call MPI_INIT (mpierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, size, mpierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD, rank, mpierr)
      tag = 100

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Lee el archivo VWM
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      if (rank .eq. 0) then
        write(*,*)"Starting vwm2rtflu_mpi ..."
        write(*,*)" Valid BC:"
        write(*,*)"  0: Neumman"
        write(*,*)"  1: U·n = 0 (Dirichlet con U como variable)"
        write(*,*)"  2: P = 0 (Dirichlet con P como variable)"
        write(*,*)" >2: Robin o Contacto"
        write(*,*)
    
c
c	abre archivo VWM
c
!      WRITE (6,*) "Archivo vwm: "
!      READ  (5,*) FILEVWM
      FILEVWM='fluido.vwm'
      OPEN  (1,FILE=FILEVWM,STATUS='OLD',ERR=210)
      write(*,*)"Abriendo fluido.vwm ..."
      
c
c     Archivo de los lados del liquido..
c
!      WRITE (6,*) "Archivo sur: "
!      READ  (5,*) FILESUR
      FILESUR='fluido.sur'
      open  (2,FILE=FILESUR,STATUS='OLD',ERR=220) ! malla de los lados (nnf)
      write(*,*)"Abriendo fluido.sur ..."
c
c	abre archivo de datos
c
!      WRITE (6,*) "Archivo de datos del preprocesador: "
!      READ  (5,*) FILEBO
      FILEBO='bc-fl.dat'
      OPEN  (3,FILE=FILEBO,STATUS='OLD',ERR=230)
      write(*,*)"Abriendo bc-fl.dat ..."

c
c	Scanea el nombre del archivo de salida
c
!      WRITE (6,*) "Archivo salida: "
!      READ  (5,*) FILEOUT
      FILEOUT='fluido.rt'
      open(12,file=FILEOUT,ERR=240)
      write(*,*)"Abriendo fluido.rt ..."

c
c     abro el archivo con las malla del liquido
c

      call searstr(1,'ELEMENT_GROUPS')
      read(1,*) el_ngrs
      nelf=0
      do 1000 i=1,el_ngrs
            read(1,*) inutil ,el_gr(i)
            nelf=nelf+el_gr(i)
 1000 continue

c
c     Leo la matriz de conectividad de los vertices debe ser <NONE>
c
      call searstr(1,'INCIDENCE')
      read(1,*)
      do 1010 i=1,nelf
         read(1,*) mmf(1,i),mmf(2,i),mmf(3,i),mmf(4,i)
 1010 continue


c
c     lee las coordenadas
c
  
      call searstr(1,'COORDINATES')
      read(1,*) nverf
      do 1020 i=1,nverf
         read(1,*) inutil,zf(1,i),zf(2,i),zf(3,i)
 1020 continue


c
c     abro el archivo con las malla del borde del liquido
c
      call searstr(2,'ELEMENT_GROUPS')
      read(2,*)sur_ngrs
      nelsur=0
      do 1030 i0=1,sur_ngrs
            read(2,*) inutil, sur_gr(i0)
            nelsur=nelsur+sur_gr(i0)
 1030 continue  

      call searstr(3,'BC_ELEMENT_GROUPS')
      do 1040 i=1,sur_ngrs
            read(3,*,end=200) NG, NR(i) ! nro de referencia y grupo
 1040 continue   

c
c     Leo la matriz de conectividad del borde
c
      call searstr(2,'INCIDENCE')
      read(2,*) dummy
      iel=1
      do 1060 i=1,sur_ngrs
            do 1050 j=1,sur_gr(i)
                  read(2,*) mmsur(1,iel),mmsur(2,iel),mmsur(3,iel)
                  mmsur(4,iel)=NR(i)
                  iel=iel+1
 1050       continue
 1060 continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c		Inicializa las matrices 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     pone a cero la matriz de caras
      do 2020 i=1, nelf
         do 2010 j=1,4
            nrcf(j,i)=0
            nnf(j,i)=0
            je(4*(i-1)+j)=mmf(j,i)
 2010    continue
 2020 continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c CONTAR LOS LADOS DE LOS TRIANGULOS Y VER LA CONECTIVIDAD!!!! 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      TI=dtime(TA)
      write(*,*)"Cálculo de matriz NNF: Starting"

      write(*,*)"trasim..."
      CALL trasim(4,je,iet,jet,nverf,nelf)
 
      endif !rank=0

      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nelf,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(nverf,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(je,nelf*4,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(iet,nverf+1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(jet,nelf*4,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)

      elemperproc = nelf / size 
      iresto = nelf - size * elemperproc
      
      if(iresto.ne.0)then
        if(rank.lt.iresto)then 
          nel_i = 1 + (elemperproc+1)*rank  
          nel_f = (elemperproc+1)*(RANK+1)
        else  
          nel_i = 1 + elemperproc*rank + iresto
          nel_f = elemperproc*(RANK+1) + iresto
        endif
      else
        nel_i = 1 + elemperproc*rank  
        nel_f = elemperproc * (RANK+1)
      endif 
      
      if(rank.eq.0) write(*,*)"conel4..."
      CALL conel4_mpi(je,iet,jet,jee,nel_i,nel_f)
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)

      do 60 i=0,size-1
        n_i=(nel_i-1)*4+1
        n_e=(nel_f-nel_i+1)*4
        CALL MPI_Bcast(n_i,1,MPI_INTEGER,i,MPI_COMM_WORLD,IERR)
        CALL MPI_Bcast(n_e,1,MPI_INTEGER,i,MPI_COMM_WORLD,IERR)
        CALL MPI_Bcast(jee(n_i),n_e,MPI_INTEGER,i,MPI_COMM_WORLD,IERR)
!        write(*,*) n_i,n_e
        if(n_e.gt.0) then
          CALL MPI_Bcast(jee(n_i),n_e,MPI_INTEGER,i,MPI_COMM_WORLD,IERR)
        endif
        CALL MPI_Barrier(MPI_COMM_WORLD,IERR)
 60   continue

      if(rank.eq.0)then

      write(*,*)"armando NNF..."
      cont=0
      do 2200 i=1,nelf
         do 2210 j=1,4
            vecino=jee(4*(i-1)+j)
            if(vecino.eq.0)then
               cont=cont+1
               nnf(j,i)=cont
               goto 2600
            else
               if(vecino.gt.i)then
                  cont=cont+1
                  nnf(j,i)=cont
                  do 2220 k=1,4
                     if(jee(4*(vecino-1)+k).eq.i)then
                        nnf(k,vecino)=cont
                        goto 2600
                     endif
 2220             continue
               endif
            endif
 2600       continue
 2210    continue
 2200 continue
      nnodf=cont
      
      write(*,*)"NFACES = ",nnodf
      write(*,*)"Cálculo de matriz NNF: Complete"

      cn(1,1)=2; cn(1,2)=3; cn(1,3)=4
      cn(2,1)=1; cn(2,2)=3; cn(2,3)=4
      cn(3,1)=1; cn(3,2)=2; cn(3,3)=4
      cn(4,1)=1; cn(4,2)=2; cn(4,3)=3

      if(nnodf.le.maxface)then
        do 2300 i=1,nnodf
           flag(i)=0
 2300   continue
        do 2400 i=1,nelf
           do 2410 j=1,4
              if(flag(nnf(j,i)).eq.0)then
                 flag(nnf(j,i))=1
                 do 2420 k=1,3
                    zbar(j,i) = (zf(k,mmf(cn(j,1),i)) + 
     &                 zf(k,mmf(cn(j,2),i)) + zf(k,mmf(cn(j,3),i)))/3
 2420            continue
              endif
 2410      continue
 2400   continue
      else
        write(*,*)"WARNING: NFACES > ",maxface
        write(*,*)"WARNING: 'fluido4.vwm' no será creado."
      endif

      TI=dtime(TA)
      write(*,*)" ", TI, "sec."

c
c	Calculo de la matriz  nrcf
c 
      TI=dtime(TA)
      write(*,*)"Calculo de matriz NRCF"

      do 3010 i=1,nelsur
         ielem=0
         N1=mmsur(1,i)
         N2=mmsur(2,i)
         N3=mmsur(3,i)
         do 3020 j=iet(N1),iet(N1+1)-1
            iev = jet(j)
            in = 4*iev-3
            do 3030 k=in,in+3
               if (je(k).eq.N2) then
                  do 3040 l=in,in+3
                     if (je(l).eq.N3) then
                        ielem=iev
                     endif
 3040             continue
               endif
 3030       continue
 3020    continue
         if (ielem.ne.0) then
            do 3100 j=1,4
               if( mmf(j,ielem).ne.N1 .and. mmf(j,ielem).ne.N2 
     &              .and. mmf(j,ielem).ne.N3 )then
                  nrcf(j,ielem)=mmsur(4,i)
                  goto 3600
               endif
 3100       continue
         endif
 3600    continue
 3010 continue

      TI=dtime(TA)
      write(*,*)" ", TI, "sec."

c
c     Archivo de salida de datos
c
c
c     Produce el archivo de salida de datos para el fluido
c
      write(12,'(a)') '*COORDINATES' ! Coordenadas
      write(12,*) nverf
      do 4000 i=1,nverf
         write(12,*) i,(zf(j,i),j=1,3)
 4000 continue
      write(12,*)

      write(12,'(/,a)') '*ELEMENT_GROUPS'
      write(12,*) el_ngrs
      do 4010 i=1,el_ngrs      
         write(12,123) i,el_gr(i),'Tetra4'
 4010 continue
      write(12,*)
 123  format(1x,i3,1x,i9,1x, a6)
      write(12,'(/,a)') '*INCIDENCE' ! Coordenadas
      write(12,'(a)') '<NONE>'
      do 4020 i=1,nelf
         write(12,*) (mmf(j,i),j=1,4)
 4020 continue
      write(12,*)

      write(12,'(/,a)') '*NFACES'
      write(12,*) nnodf
      write(12,*)

      write(12,'(/,a)') '*FACE_INCIDENCE' ! Incidencias de lados
      do 4030 i=1,nelf
         write(12,*) (nnf(j,i),j=1,4)
 4030 continue
      write(12,*)

      write(12,'(/,a)') '*NRC'    
      do 4040 i=1,nelf
         write(12,*) (nrcf(j,i),j=1,4)
 4040 continue
      write(12,*)
         
      write(12,'(/,a)') '*NSD'   !Numero de Subdominios
      write(12,*)  (' 1 ',i=1,nelf)
      write(12,*)


      write(12,'(/,a)') '*END' ! Final



c
c     Cierra los archivos de datos
c

      close(1)
      close(2)
      close(3)
      close(12)
 
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                    Creo el archivo fluido4.vwm
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(nnodf.le.maxface)then
c
c	abre archivo fluido4.vwm (para renumeracion de los nodos)
c
      OPEN  (4,FILE='fluido4.vwm',STATUS='UNKNOWN',ERR=250)

      write(4,'(A)') '*COORDINATES'
      write(4,*) nnodf
      do 6000 i=1,nnodf
            write(4,*) i,(zbar(j,i),j=1,3)
 6000 continue
      write(12,*)
      
      write(4,'(/,A)') '*ELEMENT_GROUPS'
      write(4,'(A)') ' 1'
      write(4,'(A,I7,1x,A6)') ' 1', nelf, 'Tetra4'
      write(12,*)

      write(4,'(/,A)') '*INCIDENCES'
      write(4,'(A)') '<NONE>'
      do 6100 i=1,nelf
         write(4,*) (nnf(j,i),j=1,4)
 6100 continue
      write(12,*)

      write(4,'(/,A)') '*END' ! Final
      close(4)

      endif

      endif !rank=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call MPI_FINALIZE(mpierr)
      return

 210  write (6,*) "Error al abrir el archivo: 'fluido.vwm'"
      STOP
 220  write (6,*) "Error al abrir el archivo: 'fluido.sur'"
      STOP
 230  write (6,*) "Error al abrir el archivo: 'bc-fl.dat'"
      STOP
 240  write (6,*) "Error al abrir el archivo: 'fluido.rt'"
      STOP
 250  write (6,*) "Error al abrir el archivo: 'fluido4.vwm'"
      STOP
 200  write(*,*) "Error al leer los limites de las referencias"
      STOP
      end







