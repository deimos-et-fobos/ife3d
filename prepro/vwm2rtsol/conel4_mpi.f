C     -----------------------------------------------------------------
      SUBROUTINE CONEL4_MPI (JE,IET,JET,JEE,NEL_I,NEL_F)
C     -----------------------------------------------------------------
C
      implicit double precision (a-h,o-z)
      DIMENSION JE(*), IET(*), JET(*), JEE(*)
C
      nel=nel_f-nel_i+1
      write(*,98) 'NELS per proc = ',nel
      
 98   format(A17,I8)
      iaux=int(nel/100)
      if(iaux-nel*100.ne.0) iaux=iaux+1

      DO 90 I = nel_i,nel_f
         IJE = 1 + (I-1)*4 

         N1 = JE(IJE)
         N2 = JE(IJE+1)
         N3 = JE(IJE+2)
         N4 = JE(IJE+3)
         IEV1 = 0
         IEV2 = 0
         IEV3 = 0
         IEV4 = 0
C
         DO 10 J = IET(N1),IET(N1+1)-1
            IEV = JET(J)
            IF (IEV.NE.I) THEN
               IN = 4*IEV-3
               DO 11 K=IN,IN+3
                  IF (JE(K).EQ.N2) THEN
                     DO 12 L = IN,IN+3
                        IF (JE(L).EQ.N3) THEN
                           IEV4 = IEV
                        ELSEIF (JE(L).EQ.N4) THEN
                           IEV3 = IEV
                        ENDIF
 12                  CONTINUE
                  ELSEIF (JE(K).EQ.N3) THEN
                     DO 13 L = IN,IN+3
                        IF (JE(L).EQ.N4) IEV2 = IEV
 13                  CONTINUE
                  END IF
 11            CONTINUE
            END IF
 10      CONTINUE

         DO 20 J = IET(N2),IET(N2+1)-1
            IEV = JET(J)
            IF (IEV.NE.I) THEN
               IN = 4*IEV-3
               DO 21 K=IN,IN+3
                  IF (JE(K).EQ.N3) THEN
                     DO 22 L = IN,IN+3
                        IF (JE(L).EQ.N4) IEV1 = IEV
 22                  CONTINUE
                  END IF
 21            CONTINUE
            END IF
 20      CONTINUE

         JEE(IJE)   = IEV1
         JEE(IJE+1) = IEV2
         JEE(IJE+2) = IEV3
         JEE(IJE+3) = IEV4
         
         k=I-nel_i    
         if(mod(k,iaux*10).eq.0) write(*,99) k/iaux, '%'
 99      format(2x,I2,A2)      

   90 CONTINUE

      write(*,*) '100 %'

      RETURN
      END
