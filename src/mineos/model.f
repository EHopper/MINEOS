      SUBROUTINE MODEL(IIN,IOUT)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c     read in earth model
c     write earth model to asc file
c     calc normalized elastic moduli, etc.
c     calc spline interpolation coefficients
c     calc gravity
c
c     calls: DRSPLN, GRAV
c
c     data block changed to initialization loops because common blocks must use "block data" routines, 
c     and they named to make use of these routines.  This would require renaming the unnamed common block
c     in all routines.  Init loop is an easier change, and should only be called once.
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      INTEGER*4 ITITLE(20)
      REAL*8 LCON(nknot),NCON(nknot),LSPL(3,nknot),NSPL(3,nknot)
      real*4 tdum(nknot,9)
      real*8 R(nknot),FMU(nknot),FLAM(nknot),QSHEAR(nknot),
     &       QKAPPA(nknot),XA2(nknot),XLAM(nknot),RHO(nknot),
     &       QRO(3,nknot),G(nknot),QG(3,nknot),FCON(nknot),
     &       FSPL(3,nknot),CCON(nknot),CSPL(3,nknot),ACON(nknot),
     &       ASPL(3,nknot)
      real*8 VPV(nknot),VPH(nknot),VSV(nknot),VSH(nknot),
     &            ETA(nknot),WRK(nknot10)
      COMMON R,FMU,FLAM,QSHEAR,
     &       QKAPPA,XA2,XLAM,RHO,
     &       QRO,G,QG,FCON,
     &       FSPL,LCON,LSPL,NCON,
     &       NSPL,CCON,CSPL,ACON,
     &       ASPL
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/EIFX/VPV,VPH,VSV,VSH,
     &            ETA,WRK0
     +  
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      common /head/ tdum
      DATA BIGG,TAU/6.6723D-11,1.D3/,RHOBAR/5515.D0/
c      DATA R/nknot*0.D0/,FMU/nknot*0.D0/,FLAM/nknot*0.D0/,
c     &     QSHEAR/nknot*0.D0/,QKAPPA/nknot*0.D0/,XA2/nknot*0.D0/,
c     &     XLAM/nknot*0.D0/,RHO/nknot*0.D0/,QRO/nknot3*0.D0/,
c     &     G/nknot*0.D0/,QG/nknot3*0.D0/,FCON/nknot*0.D0/,
c     &     FSPL/nknot3*0.D0/,LCON/nknot*0.D0/,LSPL/nknot3*0.D0/,
c     &     NCON/nknot*0.D0/,NSPL/nknot3*0.D0/,CCON/nknot*0.D0/,
c     &     CSPL/nknot3*0.D0/,ACON/nknot*0.D0/,ASPL/nknot3*0.D0/,
c     &     tdum/nknot9*0.0/
      PI=3.14159265358979D0
      do i=1,nknot
        r(i)=0.D0
        fmu(i)=0.D0
        flam(i)=0.D0
        qshear(i)=0.D0
        qkappa(i)=0.D0
        xa2(i)=0.D0
        xlam(i)=0.D0
        rho(i)=0.D0
        g(i)=0.D0
        fcon(i)=0.D0
        lcon(i)=0.D0
        ncon(i)=0.D0
        ccon(i)=0.D0
        acon(i)=0.D0
        do j=1,3
          qro(j,i)=0.D0
          qg(j,i)=0.D0
          fspl(j,i)=0.D0
          lspl(j,i)=0.D0
          nspl(j,i)=0.D0
          cspl(j,i)=0.D0
          aspl(j,i)=0.D0
        end do
        do j=1,9
          tdum(i,j)=0.D0
        end do
      end do

c
      READ(IIN,100) (ITITLE(I),I=1,20)
  100 FORMAT(20A4)
      READ(IIN,*) IFANIS,TREF,IFDECK
      IF(IFDECK.EQ.0) GO TO 1000
C*** CARD DECK MODEL ***
      READ(IIN,*) N,NIC,NOC
c     print*, n, nic, noc
      if (n .gt. nknot) then
        print*,' max number of knots is',nknot
        stop
      endif
      READ(IIN,105) (R(I),RHO(I),VPV(I),VSV(I),
     +     QKAPPA(I),QSHEAR(I),VPH(I),VSH(I),ETA(I),I=1,N)
c 105 FORMAT(9f10.0)
  105 format(F8.0,3F9.2,2F9.1,2F9.2,F9.5)
      GO TO 2000
C*** POLYNOMIAL MODEL ***
 1000 READ(IIN,*) NREG,NIC,NOC,RX
      RX=RX*TAU
      N=0
      KNT=0
      JJ=5
      IF(IFANIS.NE.0) JJ=8
      DO 10 NN=1,NREG
      READ(IIN,*) NLAY,R1,R2
      R1=R1*TAU
      R2=R2*TAU
      DR=(R2-R1)/FLOAT(NLAY-1)
      DO 15 I=1,NLAY
      N=N+1
   15 R(N)=R1+DR*FLOAT(I-1)
      DO 20 J=1,JJ
      READ(IIN,110) (WRK(I),I=1,5)
  110 FORMAT(5F9.5)
      DO 20 I=1,NLAY
      IND=KNT+I
      RT=R(IND)/RX
      VAL=WRK(1)+RT*(WRK(2)+RT*(WRK(3)+RT*(WRK(4)+RT*WRK(5))))
      IF(J.EQ.1) RHO(IND)=VAL*TAU
      IF(J.EQ.2) VPV(IND)=VAL*TAU
      IF(J.EQ.3) VSV(IND)=VAL*TAU
      IF(J.EQ.4) QKAPPA(IND)=VAL
      IF(J.EQ.5) QSHEAR(IND)=VAL
      IF(IFANIS.EQ.0) GOTO 20
      IF(J.EQ.6) VPH(IND)=VAL*TAU
      IF(J.EQ.7) VSH(IND)=VAL*TAU
      IF(J.EQ.8) ETA(IND)=VAL
   20 CONTINUE
   10 KNT=KNT+NLAY
 2000 IF(IFANIS.NE.0) GO TO 3000
      DO 25 I=1,N
      VPH(I)=VPV(I)
      VSH(I)=VSV(I)
   25 ETA(I)=1.D0
C*** WRITE OUT MODEL ***
 3000 WRITE(IOUT,900) (ITITLE(K),K=1,20),TREF
  900 FORMAT(1X,20A4,' REF PER =',F6.1,' SECS',///,2X,'LEVEL',
     1 4X,'RADIUS',8X,'RHO',9X,'VPV',9X,'VPH',9X,'VSV',
     2 9X,'VSH',9X,'ETA',9X,'QMU ',8X,'QKAP',/)
      WRITE(IOUT,905) (I,R(I),RHO(I),VPV(I),VPH(I),VSV(I),VSH(I),
     1 ETA(I),QSHEAR(I),QKAPPA(I),I=1,N)
c
      do i = 1, n
        tdum(i,1) = r(i)
        tdum(i,2) = rho(i)
        tdum(i,3) = vpv(i)
        tdum(i,4) = vsv(i)
        tdum(i,5) = qshear(i)
        tdum(i,6) = qkappa(i)
        tdum(i,7) = vph(i)
        tdum(i,8) = vsh(i)
        tdum(i,9) = eta(i)
      end do
  905 FORMAT(3X,I3,F12.1,5F12.2,F12.5,2F12.2)
C*** NORMALISE AND SPLINE ***
      RN=R(N)
      GN=PI*BIGG*RHOBAR*RN
      VN2=GN*RN
      VN=DSQRT(VN2)
      WN=VN/RN
      DO 44 I=1,N
      R(I)=R(I)/RN
      IF(I.GT.1.AND.DABS(R(I)-R(I-1)).LT.1.D-7) R(I)=R(I-1)
44    CONTINUE
      DO 45 I=1,N
      IF(QSHEAR(I).GT.0.D0) QSHEAR(I)=1.D0/QSHEAR(I)
      IF(QKAPPA(I).GT.0.D0) QKAPPA(I)=1.D0/QKAPPA(I)
      RHO(I)=RHO(I)/RHOBAR
      ACON(I)=RHO(I)*VPH(I)*VPH(I)/VN2
      CCON(I)=RHO(I)*VPV(I)*VPV(I)/VN2
      LCON(I)=RHO(I)*VSV(I)*VSV(I)/VN2
      NCON(I)=RHO(I)*VSH(I)*VSH(I)/VN2
      FCON(I)=ETA(I)*(ACON(I)-2.D0*LCON(I))
      FMU(I)=(ACON(I)+CCON(I)-2.D0*FCON(I)+5.D0*NCON(I)+
     1 6.D0*LCON(I))/15.D0
      FLAM(I)=(4.D0*(ACON(I)+FCON(I)-NCON(I))+CCON(I))/9.D0
     +    -2.D0*FMU(I)/3.D0
      RAT=4.D0*FMU(I)/(3.D0*(FLAM(I)+2.D0*FMU(I)))
      XLAM(I)=((1.D0-RAT)*QKAPPA(I)-.5D0*RAT*QSHEAR(I))/(1.D0-1.5D0*RAT)
   45 XA2(I)=(1.D0-RAT)*QKAPPA(I)+RAT*QSHEAR(I)
c
      CALL DRSPLN(1,N,R,RHO,QRO,WRK)
c
C*** COMPUTE G *****
c
      CALL GRAV(G,RHO,QRO,R,N)
c
      CALL DRSPLN(1,N,R,G,QG,WRK)
      CALL DRSPLN(1,N,R,FCON,FSPL,WRK)
      CALL DRSPLN(1,N,R,LCON,LSPL,WRK)
c
      IF(IFANIS.EQ.0) GOTO 60
c
      CALL DRSPLN(1,N,R,ACON,ASPL,WRK)
      CALL DRSPLN(1,N,R,CCON,CSPL,WRK)
      CALL DRSPLN(1,N,R,NCON,NSPL,WRK)
c
   60 NSL=N
      IF(VSV(NSL).GT.0.D0) GO TO 70
   65 NSL=NSL-1
      IF(VSV(NSL).LE.0.D0) GO TO 65
   70 NICP1=NIC+1
      NOCP1=NOC+1
      NSLP1=NSL+1
      TREF=0.5D0*TREF/PI
c
      RETURN
      END
