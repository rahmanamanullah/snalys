C     ****i* contourfinder/SVDCMP *
C
C   NAME
C     SVDCMP -- computes the Singular Value Decomposition (SVD).
C
C   DESCRIPTION
C     The subroutine computes the Singular Value Decomposition (SVD),
C     A = U*W*Vt for a given matrix A, and can be used to solve linear
C     equation systems if the matrix is singular.
C
C   AUTHOR
C     William H. Press, Saul A. Teukolsky, William T. Vetterling and
C     Brian P. Flannery
C
C   USAGE
C     SVDCMP( A, M, N, MP, NP, W, V )
C
C   INPUTS
C     A  - a REAL matrix with the physical dimension A(MP,NP).
C     M  - the number of rows used in A (must be less or equal to MP).
C     N  - the number of columns used in A (must be less or equal to NP).
C     MP - the number of rows in A.
C     NP - the number of columns in A.
C     W  - a REAL vector with the minimum length N.
C     V  - a REAL matrix with the dimensions V(N,N).
C
C   SIDE EFFECTS
C     A contains the orthonormal U matrix on return. W and V contains the
C     corresponding vector and matrix in the method.
C
C   SEE ALSO
C     Numerical Recepies - The Art of Scientific Computing, Second Edition,
C     by William H. Press, Saul A. Teukolsky, William T. Vetterling and
C     Brian P. Flannery, section 2.6 Singular Value Decomposition.
C
C   USES
C     PYTHAG()
C
C   USED BY
C     NEXTNETPOINT()
C
C     ***
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      INTEGER M,MP,N,NP,NMAX
      REAL A(MP,NP),V(NP,NP),W(NP)
      PARAMETER (NMAX=500)
CU    USES pythag
      INTEGER I,ITS,J,JJ,K,L,NM
      REAL ANORM,C,F,G,H,S,SCALE,X,Y,Z,RV1(NMAX),PYTHAG
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF(I.LE.M)THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF(SCALE.NE.0.0)THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            DO 15 J=L,N
              S=0.0
              DO 13 K=I,M
                S=S+A(K,I)*A(K,J)
13            CONTINUE
              F=S/H
              DO 14 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
14            CONTINUE
15          CONTINUE
            DO 16 K=I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF((I.LE.M).AND.(I.NE.N))THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF(SCALE.NE.0.0)THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            DO 23 J=L,M
              S=0.0
              DO 21 K=L,N
                S=S+A(J,K)*A(I,K)
21            CONTINUE
              DO 22 K=L,N
                A(J,K)=A(J,K)+S*RV1(K)
22            CONTINUE
23          CONTINUE
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF(I.LT.N)THEN
          IF(G.NE.0.0)THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=MIN(M,N),1,-1
        L=I+1
        G=W(I)
        DO 33 J=L,N
          A(I,J)=0.0
33      CONTINUE
        IF(G.NE.0.0)THEN
          G=1.0/G
          DO 36 J=L,N
            S=0.0
            DO 34 k=l,m
              S=S+A(K,I)*A(K,J)
34          CONTINUE
            F=(S/A(I,I))*G
            DO 35 K=I,M
              A(K,J)=A(K,J)+F*A(K,I)
35          CONTINUE
36        CONTINUE
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF((ABS(RV1(L))+ANORM).EQ.ANORM)  GOTO 2
            IF((ABS(W(NM))+ANORM).EQ.ANORM)  GOTO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            RV1(I)=C*RV1(I)
            IF((ABS(F)+ANORM).EQ.ANORM) GOTO 2
            G=W(I)
            H=PYTHAG(F,G)
            W(I)=H
            H=1.0/H
            C= (G*H)
            S=-(F*H)
            DO 42 J=1,M
              Y=A(J,NM)
              Z=A(J,I)
              A(J,NM)=(Y*C)+(Z*S)
              A(J,I)=-(Y*S)+(Z*C)
42          CONTINUE
43        CONTINUE
2         Z=W(K)
          IF(L.EQ.K)THEN
            IF(Z.LT.0.0)THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GOTO 3
          ENDIF
          IF(ITS.EQ.30) PAUSE 'no convergence in svdcmp'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=PYTHAG(F,1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=PYTHAG(F,H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 JJ=1,N
              X=V(JJ,J)
              Z=V(JJ,I)
              V(JJ,J)= (X*C)+(Z*S)
              V(JJ,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=PYTHAG(F,H)
            W(J)=Z
            IF(Z.NE.0.0)THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 JJ=1,M
              Y=A(JJ,J)
              Z=A(JJ,I)
              A(JJ,J)= (Y*C)+(Z*S)
              A(JJ,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
C     *** END OF SUBROUTINE SVDCMP() ***

C     END OF FILE
