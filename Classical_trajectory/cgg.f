      IMPLICIT REAL*8(A-H,O-Z)
C--------------------------------------------
      PARAMETER (MATWO=27)
      PARAMETER (MATHR=169)
      COMMON/COEFF/COFTWO(MATWO),COFTHR(MATHR)
C--------------------------------------------
      COMMON/AMASS/AM1,AM2,AM3
      COMMON/MASS/EM(3),XM,EMU,DM(3),EPS2,EPS3
      COMMON/REA/ICH 
      COMMON/VPOT/EN,HH0,EEK,EC1,EC2,AMYS
      COMMON/AIN/AA1,AA2,AA3
      COMMON/SPEC/BE,CE,AE,REE,AV1,DDD,BETA  
      COMMON/ENERGY/EKINA,EKINR  
      COMMON/PPHI/ETOTA,VPOTA
      common/arand/ix,iy,iz
      COMMON/ENER/EVEC,HSTRE
      COMMON/FACT/SWT1,SWT2
      COMMON/RJL/RJA,RLA
      DIMENSION Y(12),D(12),RAN(9)
C     DIMENSION E(12),W1(12),W2(12),W3(12),W4(12)  
      DIMENSION RE(3),SUM(3),OM(3),DD(3),BET(3)
      DIMENSION G(12),DERY(12),YP(12),DAR(12,15),DP(12,15)
      DIMENSION RA(3),GEM(3,10,50),HYPJ(3)
      DIMENSION SUMMJ(3,10),SUMMV(3,50)
      DIMENSION PRB1(3,10,50),PRB2(3,10,50) 
      DATA PI/3.14159265359D0/
      EXTERNAL DERCL 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      OPEN(1,FILE='TWOBODY_COFF.dat',STATUS='OLD')
      OPEN(2,FILE='THREEBODY_COFF.dat',STATUS='OLD')
*
      DO I=1,MATWO
      READ(1,*)COFTWO(I)
      ENDDO
*
      DO I=1,MATHR
      READ(2,*)COFTHR(I)
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      READ(5,100) RE(1),RE(2),RE(3) ! Three equlibium diatomic bond legths 
      READ(5,100) AM1,AM2,AM3       ! Three masses of atoms
      READ(5,100) DD(1),BET(1),DD(2),BET(2),DD(3),BET(3) ! Three pairs of Morse Parametres for the diatoms 
      READ(5,200) N1,IBG,IX
      READ(5,*) R,RMAX,BMAX,EDEF ! EDEF is the translational energy [1.0 ev = 8065.5 cm^(-1)]
      READ(5,201) PRMT
      READ(5,203) NTR,NMAX,NBOX,JMAX,JBOX
      READ(5,200) IOUT
      READ(5,200) IANA
      READ(5,*) SWT1,SWT2
  100 FORMAT(6F6.1)  
  200 FORMAT(2I3,I10)
  201 FORMAT(E7.1)
  203 FORMAT(I6,4I3)
 
C IF IANA=1 ANALYSES OF DATA ON UNIT 9.
C IF IOUT=1 DETAILED OUTPUT ALONG TRAJECTORY 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      AMR=AM2*AM3/(AM2+AM3) ! Diatomic (2-3) reduced mass 
      AMY=AM1*(AM2+AM3)/(AM1+AM2+AM3) !Triatomic center of mass 
      AM=AM1+AM2+AM3 ! Total Mass
      AMYS=DSQRT(AM1*AM2*AM3/AM) ! Three body reduced mass 
      DM(1)=DSQRT(AM1*(AM2+AM3)/AMYS/AM)
      DM(2)=DSQRT(AM2*(AM1+AM3)/AMYS/AM)
      DM(3)=DSQRT(AM3*(AM1+AM2)/AMYS/AM)
      EPS2=2.0*DATAN(AM3/AMYS)
      EPS3=2.0*DATAN(AM2/AMYS)
CCC   1 cm^(-1)= 0.0001196265728 epsilon
CCC   EK is the translational energy in epsilon unit  
      EK=0.0001196265728D0*EDEF 
      HSTRE=0.06350781278D0  ! hcross (amu.ang**2/tau) 
      VX=SQRT(2.0D0*EK/AMYS) ! velocity 
c     ALL=AMYS*VX*BMAX/HSTRE-0.5
c     WRITE(6,*)'ALL=',all 
c     WRITE(6,*)'BMAX=',BMAX
      BMAX=HSTRE*10.5/(AMYS*VX)
c     WRITE(6,*)'BMAX=',BMAX
      IF (IOUT.EQ.1) WRITE(6,300) EDEF,EK
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Initialization of Coordinates and Momenta
C JCP, 63, 2214 (1975), R. N. Porter, L. M. Raff, W. H. Miller
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RE1=RE(1)
      BETA1=BET(1)
      D1=DD(1) 
      RHOFI=R+RE1+0.3
      AI1=AMR*RE1*RE1
      RN1=1.*N1
      RCCC=1./BETA1/RE1
      RIM=AMR*RE1*RE1
      AAVR=(1.-3.*RCCC*(1.-RCCC))/2./RIM
      BBVR=2.*RCCC*(1.-1.5*RCCC)/RIM
      CCVR=RCCC*(1.-3.*RCCC)/2./RIM
      RHST=1./HSTRE  
      HST2=0.5*HSTRE**2
      BE1=HST2/AMR/RE1/RE1 
      DE1P=HST2**2/AMR**2/BETA1**2/RE1**6/D1 
      XE1=HSTRE*BETA1/DSQRT(8.*AMR*D1) 
      HME1=HSTRE*BETA1*DSQRT(2.*D1/AMR)
      ALE1=1.5*HST2*HME1*(1.-1./BETA1/RE1)/AMR/BETA1/RE1**3/D1 
      EVR1I=(1.*IBG+0.5)**2*(BE1-DE1P*(1.*IBG+0.5)**2-ALE1*(RN1+0.5))
     1     +HME1*(RN1+0.5)*(1.-XE1*(RN1+0.5))
      AK1=1./XE1
      EV1=D1-D1*(1.-(2.*RN1+1.)/AK1)**2
      RJ1=1.0*IBG+0.5
      EROTA=RJ1**2*(BE1-(RN1+0.5)*ALE1)-RJ1**4*DE1P
      EVR1I=RJ1**2*(BE1-DE1P*RJ1**2-ALE1*(RN1+0.5)) 
     1     +HME1*(RN1+0.5)*(1.-XE1*(RN1+0.5)) ! E_vJ in Eq.48, W. H. Miller, JCP, 63, 2214 (1975)
      ALAM1=HSTRE*RJ1
      AL21=ALAM1**2  
      AVR1=EVR1I-D1-AAVR*AL21 ! a in Eq.6, W. H. Miller, JCP, 63, 2214 (1975)
      BVR1=2.*D1-BBVR*AL21    ! b in Eq.6, W. H. Miller, JCP, 63, 2214 (1975) 
      CVR1=-D1+CCVR*AL21      ! c in Eq.6, W. H. Miller, JCP, 63, 2214 C(1975)
      ALA1=2.*DSQRT(AVR1*CVR1)/BVR1
      ALA12=-2.*AVR1/BVR1  
      IF (ALA1.GT.1.0) ALA1=1.0D0
      IXG=IX
      ALA1S=DSQRT(DABS(1.-ALA1**2))
  220 CONTINUE 
      EK=0.0001196265728D0*EDEF  
  123 FORMAT(1X,F12.8)
      DO 10 K=1,9
      YFL=random(0)
      RAN(K)=YFL
   10 CONTINUE 
      QVA1=6.283185*RAN(5) 
      QVA2=6.283185*RAN(6) 
      COG1=DCOS(QVA1)
      QI1=ALA12/(1.-ALA1S*COG1) ! Xi in Eq.36, W. H. Miller, JCP, 63, 2214 (1975)  
      R1A=RE1-DLOG(QI1)/BETA1   ! r_e in Eq.35, W. H. Miller, JCP, 63, 2214 (1975)
      PR1=DSQRT(2.*DABS(AVR1+BVR1*QI1+CVR1*QI1**2)/AMR) ! p_r in Eq.5, W. H. Miller, JCP, 63, 2214 (1975)  
      DELA1=DSQRT(-0.5*AMR/CVR1)*(BVR1*CCVR/CVR1+2.*BBVR)/BETA1*
     1(DASIN((QI1-ALA12)/QI1/ALA1S)+DASIN((1.+2.*CVR1*QI1/BVR1)
     2/ALA1S))+AMR*PR1*CCVR/CVR1/BETA1 ! Delta_l in Eq.26, W. H. Miller, JCP, 63, 2214 (1975) 
      ITER=0
  503 V1=D1*(DEXP(-BETA1*(R1A-RE1))-1.)**2 ! Morse Potential in Eq.1, W. H. Miller, JCP, 63, 2214 (1975)
      PR1=EVR1I-V1-0.5*AL21/AMR/R1A**2 
      ITER=ITER+1
      IF (PR1.GE.0.00) GO TO 501 
      R1A=R1A+(RE1-R1A)/DABS(RE1-R1A)*0.01
      IF (ITER.LT.2) GO TO 503
      IF (ITER.LT.50) GO TO 503  
  501 PR1=DSQRT(2.*DABS(PR1)/AMR)
      IF (QVA1.GT.3.141592654) PR1=-PR1
      A1=3.14159265*RAN(1) 
      B1=6.283185*RAN(2)
      C1=6.283185*RAN(3)
      EK=EK-HST2*RJ1**2/AI1
      RLM=RHST*DSQRT(2.*AMY*EK)*BMAX
      RL=RLM*RAN(4)+0.5
      E1=0.0001196265728*EDEF-EROTA
      VR=-DSQRT(2.*E1/AMY) 
      FAC=3.14159265*2.*HST2*RL*RLM/AMY/E1
      IF (IANA.EQ.1) READ(5,141) IX,FAC
      IF (IANA.EQ.1) GO TO 4000  
C     IF (JUMP.EQ.1) RL=BMAX*RHST*AMY*DABS(VR)
      BIMPA=RL/AMY/RHST/DABS(VR) 
      QL=6.283185*RAN(7)
      A3=3.14159265*RAN(8) 
      C3=6.283185*RAN(9)
      FACTM=6.28319*RAN(4)*BMAX**2
      EORB=HST2*RL**2/AMY/R/R
      EK1=E1-EORB
      SA3=DSIN(A3)
      CA3=DCOS(A3)
      SC3=DSIN(C3)
      CC3=DCOS(C3)
      CQL=DCOS(QL)
      SQL=DSIN(QL)
      B1=B1+DELA1*ALAM1
      SB1=DSIN(B1)
      CB1=DCOS(B1)
      SC1=DSIN(C1)
      CC1=DCOS(C1)
      SA1=DSIN(A1)
      CA1=DCOS(A1)
      R1=-CA3*SQL*SC3+CQL*CC3
      R2=CA3*SQL*CC3+CQL*SC3
      R3=-SA3*SQL
      Y(1)=R*R1
      Y(2)=R*R2
      Y(3)=R*R3
      FA=HSTRE*RL/R/AMY
      Y(7)=R1A*(CB1*CC1-CA1*SB1*SC1)
      Y(8)=R1A*(CB1*SC1+CA1*SB1*CC1)
      Y(9)=-R1A*SA1*SB1
C CORRECTION FOR FINITE POTENTIAL AT DISTANCE CHOSEN
      RA(2)=R1A
      RA(3)=DSQRT((AM2*Y(1)+AM3*(Y(1)+Y(7)))**2+
     1(AM2*Y(2)+AM3*(Y(2)+Y(8)))**2+(AM2*Y(3)+AM3*(Y(3)+Y(9)))**2)
     1/(AM2+AM3)
      RA(1)=DSQRT((AM3*Y(1)+AM2*(Y(1)-Y(7)))**2+
     1(AM3*Y(2)+AM2*(Y(2)-Y(8)))**2+(AM3*Y(3)+AM2*(Y(3)-Y(9)))**2)
     1/(AM2+AM3)
      DO 50 I=1,3
   50 RA(I)=RA(I)/0.5291770644
      CALL POTENTIAL(RA,EN)
C     CALL POTH3(RA,EN) 
      EN=EN*0.9648455078
      EK1=EK1-D1-EN+V1
C     WRITE(6,400) EK1,EN,V1,EORB,EROTA,D1
C     IF (EK1.LT.0.0) GO TO 222  
      VR=-DSQRT(2.*EK1/AMY)
      Y(4)=VR*R1+FA*(-CA3*CQL*SC3-SQL*CC3)
      Y(5)=VR*R2+FA*(CA3*CQL*CC3-SQL*SC3)
      Y(6)=VR*R3-FA*CQL*SA3
      FA=ALAM1/R1A/AMR
      Y(10)=FA*(-CA1*CB1*SC1-SB1*CC1)+PR1*(-CA1*SB1*SC1+CB1*CC1)
      Y(11)=FA*(CA1*CB1*CC1-SB1*SC1)+PR1*(CA1*SB1*CC1+CB1*SC1) 
      Y(12)=-FA*CB1*SA1-PR1*SA1*SB1
C     WRITE(6,300) (Y(I),I=1,12) 
      ICH=1
      EKINA=0.5*AMR*(Y(10)**2+Y(11)**2+Y(12)**2)
      EKINR=0.5*AMY*(Y(4)**2+Y(5)**2+Y(6)**2)
      VPOTA=EN 
      ETOTA=EKINA+EKINR+EN 
C     WRITE(6,2000) EKINA,EKINR,EN,ETOTA
 2000 FORMAT(' ETEST ',4F10.6)
      CALL CARDEL(Y,D)
C USED IN CONNECTION WITH D Y(12)/DT=0 
C     WRITE(6,300) (D(I),I=1,12) 
C     IF (D(1).GT.RHOFI) WRITE(6,4001) D(1)  
      IF (D(1).GT.RHOFI) STOP
 4001 FORMAT(' INITIAL RHO=',F10.5)
      DO 309 I=1,12  
  309 YP(I)=D(I)
      CALL DELCAR(YP,Y)
      EKINA=0.5*AMR*(Y(10)**2+Y(11)**2+Y(12)**2)
      EKINR=0.5*AMY*(Y(4)**2+Y(5)**2+Y(6)**2)
      RA(2)=DSQRT(Y(7)**2+Y(8)**2+Y(9)**2)
      RA(3)=DSQRT((AM2*Y(1)+AM3*(Y(1)+Y(7)))**2+
     1(AM2*Y(2)+AM3*(Y(2)+Y(8)))**2+(AM2*Y(3)+AM3*(Y(3)+Y(9)))**2)
     1/(AM2+AM3)
      RA(1)=DSQRT((AM3*Y(1)+AM2*(Y(1)-Y(7)))**2+
     1(AM3*Y(2)+AM2*(Y(2)-Y(8)))**2+(AM3*Y(3)+AM2*(Y(3)-Y(9)))**2)
     1/(AM2+AM3)
      DO 59 I=1,3
   59 RA(I)=RA(I)/0.5291770644
      CALL POTENTIAL(RA,EN)
C     CALL POTH3(RA,EN) 
      EN=EN*0.9648455078
      ETOTA=EKINA+EKINR+EN 
C     WRITE(6,2000) EKINA,EKINR,EN,ETOTA
  305 FORMAT(' INVERSE TEST')
C     WRITE(6,300) (Y(I),I=1,12) 
C START INTEGRATION OF TRAJECTORY
      IF (D(7).GT.0.0) D(7)=-D(7)
C     WRITE(6,202) EDEF,BIMPA
  202 FORMAT(' KIN. ENERGY/CM-1+IMP.PARAM. ',2E10.3)
      IFAIL=0  
      IT=2
      N=12
      N1=12
      IQ=0
      X=0.0
      H=0.10
      H0=H
      DO 18 I=1,N
   18 G(I)=PRMT
      CALL DERCL(DERY,D,X) 
C     CALL DYBDX(D,DERY,N,X)
      ET=EN+HH0+EEK+EC1+EC2+EVEC
C     WRITE(6,409) EEK,EN,HH0,EC1,EC2,ET
  409 FORMAT(' INITIAL ENERGIES ',6F10.6)
   20 CALL D02AHF(X,D,G,IT,N,H0,H,DERCL,DERY,YP,DAR,DP,N1,IQ,IFAIL)  
C  20 CALL DA01AD(D,E,W1,W2,W3,W4,N,X,H)
C CALCULATION OF ORBITAL AND ROTATIONAL ANGULAR MOMENTUM 
      AMRO=0.25*AMYS*D(1)*D(1)
      SB=DSIN(D(5))  
      CB=DCOS(D(5))  
      CG=DCOS(D(6))  
      SG=DSIN(D(6))  
      AJX=-D(10)*CG*SB
      AJY=-D(10)*SG*SB
      AJZ=-D(10)*CB  
      CP=DCOS(D(3))  
      SP=DSIN(D(3))  
      CT=DCOS(D(2))  
      ST=DSIN(D(2))  
      RJX=0.5*(AJX*(1.-CP)-AJY*CT*SP/(1.+ST))
      RJY=0.5*(AJY*(1.+CP)-AJX*CT*SP/(1.-ST))
      RLX=0.5*(AJX*(1.+CP)+AJY*CT*SP/(1.+ST))
      RLY=0.5*(AJY*(1.-CP)+AJX*CT*SP/(1.-ST))
      RJZ=-D(8)*SP+0.5*AJZ+CP*(AJZ-2.*CT*D(9))/ST  
      RLZ=D(8)*SP+0.5*AJZ-CP*(AJZ-2.*CT*D(9))/ST
      RJTX=RJX+RLX
      RJTY=RJY+RLY
      RJTZ=RJZ+RLZ
      RJTOT=DSQRT(RJTX**2+RJTY**2+RJTZ**2)
      IF (IOUT.EQ.1) WRITE(6,1009) RJX,RJY,RJZ,RLX,RLY,RLZ,RJTOT
 1009 FORMAT(' RJTOT=',7E10.3)
      ET=EN+HH0+EEK+EC1+EC2+EVEC
      IF (IOUT.EQ.1) WRITE(6,423) X,D(1),D(2),D(3),EEK,EN,HH0  
     1,EC1,EC2,EVEC,ET
C     WRITE(6,400) (D(I),I=4,12) 
      DO 444 I=1,12  
  444 YP(I)=D(I)
      IF (IOUT.EQ.1) CALL DELCAR(YP,Y) 
C     PIH=0.5*PI
C     H2=2.*H  
C     IF (DABS(D(2)-PIH).LT.H2) H=0.5*H
C     IF (DABS(D(2)-PIH).GT.H2) H=2.*H 
      IF (H.GT.H0) H=H0
  423 FORMAT(1X,11E10.3)
  400 FORMAT(1X,10E10.3)
      IF (D(1).LT.RHOFI) GO TO 20
 4000 CONTINUE 
      IF (IANA.EQ.1) READ(5,140) X,(D(I),I=1,12)
      HYPJ(1)=D(1)
      HYPJ(2)=D(2)
      HYPJ(3)=D(3)
      WRITE(6,141) IXG,FAC,BMAX,EDEF
  141 FORMAT(1X,I10,6F12.5)
      WRITE(6,140) X,(D(I),I=1,12)
  140 FORMAT(1X,5D15.9)
      CALL GETR(HYPJ,RA)
      IF (IOUT.EQ.1) WRITE(6,400) (RA(I),I=1,3)
      IF (RA(2).GT.RMAX.AND.RA(3).GT.RMAX) ICH=3
      IF (RA(1).GT.RMAX.AND.RA(3).GT.RMAX) ICH=1
      IF (RA(1).GT.RMAX.AND.RA(2).GT.RMAX) ICH=2
  500 FORMAT(' REACTIVE')  
      CALL DELCAR(D,Y)
      IF (IOUT.EQ.1) PRINT 305
      IF (IOUT.EQ.1) WRITE(6,300) (Y(I),I=1,12)
      AA1=AM2  
      AA2=AM3  
      AA3=AM1  
      IF (ICH.EQ.1) GO TO 301
      AA1=AM1  
      AA3=AM2  
      IF (ICH.EQ.2) GO TO 301
      AA2=AM2  
      AA3=AM3  
  301 CONTINUE 
C     WRITE(6,302) ICH
  302 FORMAT(' CHANNEL=',I3)
      AMM=AA1*AA2/(AA1+AA2)
      OM(ICH)=BET(ICH)*DSQRT(2.*DD(ICH)/AMM) 
      RE2=RE(ICH)*RE(ICH)  
      BE=.00201662/AMM/RE2 
      CE=4.0668E-6/(AMM*RE2*BET(ICH))**2/RE2/DD(ICH)
      AE=1.921068E-4*OM(ICH)*(1.-1./(BET(ICH)*RE(ICH)))  
     1/AMM/BET(ICH)/RE(ICH)**3/DD(ICH) 
      AV1=2.*DSQRT(2.*AMM*DD(ICH))*15.74609416/BET(ICH)  
      BETA=BET(ICH)  
      REE=RE(ICH)
      DDD=DD(ICH)
C     WRITE(6,987) BE,CE,AE,AV1,BETA,REE,DDD 
  987 FORMAT(' SPECTR. BE,CE,AE,AV1,BETA,RE,D ',7E10.3)  
      CALL ROTVIB(Y,EROT,EVIB,RJA,AV,RL)
      ETOTA=EKINA+EKINR+EN 
C     WRITE(6,2000) EKINA,EKINR,EN,ETOTA
      ITR=ITR+1
C     WRITE(6,303) ITR,EROT,EVIB,RJA,AV,RLA
  303 FORMAT(' ITR,EROT,EVIB,J,V,L ',I4,5E10.3)
      DO 401 I=1,NMAX
      VMI=1.*(I-1)-0.5*NBOX
      VMA=VMI+1.*NBOX
      IF (AV.GT.VMI.AND.AV.LE.VMA) N=I 
      DO 401 J=1,JMAX,JBOX 
      RMI=1.*(J-1)-0.5*JBOX
      RMA=RMI+1.*JBOX
      IF (RJA.GT.RMI.AND.RJA.LE.RMA) JJ=J
  401 CONTINUE 
      GEM(ICH,N,JJ)=GEM(ICH,N,JJ)+FAC  
      IF (ITR.LT.NTR) GO TO 220  
      DO 402 I=1,3
      DO 402 J=1,NMAX
      DO 402 K=1,JMAX
  402 GEM(I,J,K)=GEM(I,J,K)/NTR  
      IF (IOUT.EQ.1) WRITE(6,406) NTR  
  406 FORMAT(' NTR=',I4)
      IF (IOUT.EQ.1) PRINT 403
  403 FORMAT('   N  J  CROSS SECTIONS IN AA**2 FOR ICH=1,2,3') 
      RNTR=FLOAT(NTR)
      DO 404 I=1,NMAX
      I1=I-1
      DO 404 J=1,JMAX,JBOX 
      J1=J-1
      IF (IOUT.EQ.1) WRITE(6,405) I1,J1,(GEM(K,I,J),K=1,3)
     1,GEM(2,I,J)+GEM(3,I,J)
      SUM(1)=SUM(1)+GEM(1,I,J)
      SUM(2)=SUM(2)+GEM(2,I,J)
      SUM(3)=SUM(3)+GEM(3,I,J)
      RNIJ1=GEM(1,I,J)*RNTR
      RNIJ2=GEM(2,I,J)*RNTR
      RNIJ3=GEM(3,I,J)*RNTR
      PRB1(1,I,J)=(1.0D0/RNTR)
     1*(RNIJ1+RNIJ1*SQRT(1.0D0/RNIJ1-1.0D0/RNTR))
      PRB1(2,I,J)=(1.0D0/RNTR)
     1*(RNIJ2+RNIJ2*SQRT(1.0D0/RNIJ2-1.0D0/RNTR))
      PRB1(3,I,J)=(1.0D0/RNTR)
     1*(RNIJ3+RNIJ3*SQRT(1.0D0/RNIJ3-1.0D0/RNTR))
      PRB2(1,I,J)=(1.0D0/RNTR)
     1*(RNIJ1-RNIJ1*SQRT(1.0D0/RNIJ1-1.0D0/RNTR))
      PRB2(2,I,J)=(1.0D0/RNTR)
     1*(RNIJ2-RNIJ2*SQRT(1.0D0/RNIJ2-1.0D0/RNTR))
      PRB2(3,I,J)=(1.0D0/RNTR)
     1*(RNIJ3-RNIJ3*SQRT(1.0D0/RNIJ3-1.0D0/RNTR))
  404 CONTINUE 
CCCC
      IF (IOUT.EQ.1) THEN
        DO K=1,3
          DO I=1,NMAX
            SUMMJ(K,I)=0.0D0
            DO J=1,JMAX,JBOX
              SUMMJ(K,I)=SUMMJ(K,I)+GEM(K,I,J)
            ENDDO  
          ENDDO  
*
          DO J=1,JMAX,JBOX
            SUMMV(K,J)=0.0D0
            DO I=1,NMAX
              SUMMV(K,J)=SUMMV(K,J)+GEM(K,I,J)
            ENDDO  
          ENDDO  
        ENDDO  
*
        DO I=1,NMAX
          I1=I-1
          WRITE(10,99)I1,(SUMMJ(K,I),K=1,3),SUMMJ(2,I)+SUMMJ(3,I) 
        ENDDO  
*
        DO J=1,JMAX,JBOX
          J1=J-1
          WRITE(20,99)J1,(SUMMV(K,J),K=1,3),SUMMV(2,J)+SUMMV(3,J) 
        ENDDO  
      ENDIF 
   99  FORMAT(1X,I3,4F12.6)
c     IF (IOUT.EQ.1) WRITE(6,*) 'POSITIVE ERROR' 
c     DO 563 I=1,NMAX
c     I1=I-1
c     DO 563 J=1,JMAX,JBOX 
c     J1=J-1
c     IF (IOUT.EQ.1) WRITE(6,405) I1,J1,(PRB1(K,I,J),K=1,3)
c 563 CONTINUE 
c     IF (IOUT.EQ.1) WRITE(6,*) 'NEGATIVE ERROR' 
c     DO 564 I=1,NMAX
c     I1=I-1
c     DO 564 J=1,JMAX,JBOX 
c     J1=J-1
c     IF (IOUT.EQ.1) WRITE(6,405) I1,J1,(PRB2(K,I,J),K=1,3)
c 564 CONTINUE 
  405 FORMAT(1X,2I3,4F10.5)
      IF (IOUT.EQ.1) WRITE(6,407) (SUM(I),I=1,3),SUM(2)+SUM(3)
  407 FORMAT(' TOTAL ',4F10.5)
  300 FORMAT(1X,10E10.3)
      STOP
      END
      SUBROUTINE GETR(HYPJ,R)
C  
C  CONVERT JOHNSON'S HYPERSPHERICAL COORDINATES TO INTERNAL DISTANCE 
C  COORDINATES (R1,R2,R3)  
C  
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      DIMENSION HYPJ(3),R(3)
      COMMON/MASS/EM(3),XM,EMU,DM(3),EPS2,EPS3
      ST=DSIN(HYPJ(2))
      FAC=HYPJ(1)/DSQRT(2.D0)
      R(1)=FAC*DM(3)*DSQRT(1.D0+ST*DCOS(HYPJ(3)+EPS3))
      R(2)=FAC*DM(1)*DSQRT(1.D0+ST*DCOS(HYPJ(3)))  
      R(3)=FAC*DM(2)*DSQRT(1.D0+ST*DCOS(HYPJ(3)-EPS2))
      RETURN
C  END OF GETR 
      END
CCCC  This Routine contains the expression of derivatives of Y
CCCC  i.e. The differential equations of six coordinates and six momentums
      SUBROUTINE DERCL(DERY,Y,X) 
C     SUBROUTINE DYBDX(Y,DERY,N,X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MASS/EM(3),XM,EMU,DM(3),EPS2,EPS3
      COMMON/VPOT/EN,H0,EK,EC1,EC2,AMYS
      DIMENSION Y(12),DERY(12),YDOT(12)
      DIMENSION HYPJ(3),DV(3),R(3)
      DATA EPSH,A0/.9648455078D0,0.5291770644D0/
      COMMON/COEFF/COFTWO(27),COFTHR(169)
*
      HYPJ(1)=Y(1) ! Rho
      HYPJ(2)=Y(2) ! Theta
      HYPJ(3)=Y(3) ! Phi
      CALL GETR(HYPJ,R)
      DO I=1,3
      R(I)=R(I)/A0
      ENDDO
      CALL POTENTIAL(R,EN)
C     CALL POTH3(R,EN)  
      EN=EN*EPSH
      CALL DERIV_POT(R,DV) ! R in bohr, DV in ev/bohr
      DO I=1,3
      R(I)=R(I)*A0
      DV(I)=DV(I)*EPSH/A0 ! DV in eps/ang
      ENDDO
C     WRITE(90,90)R(1),R(2),R(3),EN
C     WRITE(91,90)R(1),R(2),R(3),DV(1),DV(2),DV(3)
C  90 FORMAT(1X,6F12.6)
CCCC  
CCCC JCP, 79, 1906 (1983), Johnson: Appendix B
CCCC  
      DERY(1)=Y(7)/AMYS ! derivative of hyperradius (rho)
      AMRO=AMYS*Y(1)*Y(1)  
      ST=DSIN(Y(2))  
      CT=DCOS(Y(2))  
      ST2=ST*ST
      H0=2.*(Y(8)**2+Y(9)**2/ST2)/AMRO 
      PARA=1.0D0+ST*DCOS(Y(3)+EPS3)
      PARB=1.0D0+ST*DCOS(Y(3))
      PARC=1.0D0+ST*DCOS(Y(3)-EPS2)
      DR1DR=R(1)/Y(1)
      DR2DR=R(2)/Y(1)
      DR3DR=R(3)/Y(1)
      FT=0.5*Y(1)/DSQRT(2.D0)
      DR1DT=FT*DM(3)*CT*DCOS(Y(3)+EPS3)/DSQRT(PARA)
      DR2DT=FT*DM(1)*CT*DCOS(Y(3))/DSQRT(PARB)
      DR3DT=FT*DM(2)*CT*DCOS(Y(3)-EPS2)/DSQRT(PARC)
      DR1DP=-FT*DM(3)*ST*DSIN(Y(3)+EPS3)/DSQRT(PARA)
      DR3DP=-FT*DM(2)*ST*DSIN(Y(3)-EPS2)/DSQRT(PARC)
      DR2DP=-FT*DM(1)*ST*DSIN(Y(3))/DSQRT(PARB)
      CT2=CT*CT
      DERY(2)=4.*Y(8)/AMRO ! derivative of theta (hyperangle), Y(8)= P_theta
      DERY(3)=(4.*Y(9)-2.*CT*Y(12))/AMRO/ST2 ! derivative of phi (hyperangle), Y(9)= P_phi, Y(12)= P_gamma
      PAR=1.-ST*DCOS(2.*Y(6))
      DERY(4)=2.*Y(10)*PAR/AMRO/CT2 ! derivative of alpha, Y(10)= P_alpha, Y(6)= gamma
      DERY(5)=0.0 ! derivative of beta 
      DERY(6)=(Y(12)-2.*CT*Y(9))/AMRO/ST2-2.*Y(12)*PAR/AMRO/CT2 ! derivative of gamma
      DERY(10)=0.0 ! derivative of P_alpha
      DERY(11)=0.0 ! derivative of P_beta
      AMROCT=AMRO*CT2
      PAR2=Y(10)**2-Y(12)**2
      DERY(9)=-DV(1)*DR1DP-DV(2)*DR2DP-DV(3)*DR3DP ! derivative of P_phi
      DERY(8)=-DV(1)*DR1DT-DV(2)*DR2DT-DV(3)*DR3DT  
     1+((4.*Y(9)**2*CT/ST2-2.*Y(12)*Y(9))/ST+Y(12)*(Y(12)-4.*CT*
     1Y(9))*CT/ST2/ST-2.*PAR2*ST*PAR/CT2/CT+PAR2*DCOS(2.*Y(6)) 
     1/CT)/AMRO                                    ! derivative of P_theta
      DERY(7)=-DV(1)*DR1DR-DV(2)*DR2DR-DV(3)*DR3DR 
     1+(Y(12)*(Y(12)-4.*CT*Y(9))/ST2+2.*PAR2*PAR/CT2)/AMRO/Y(1)
     1+2.*H0/Y(1)                                  ! derivative of P_rho
      DERY(12)=-2.*ST*DSIN(2.*Y(6))*PAR2/AMROCT    ! derivative of P_gamma

CCCC  VECTOR POTENTIAL TERM in Eq.7: JCP, 107, 6213 (1997) Adhikari & Billing
C     CALL VECP(Y,AMYS,YDOT)
C     DO 3 I=1,12
C     DERY(I)=DERY(I)+YDOT(I)
C  3  CONTINUE

      EK=0.5*Y(7)**2/AMYS  
      EC1=0.5*Y(12)*(Y(12)-4.*CT*Y(9))/AMRO/ST2
      EC2=PAR2*PAR/AMRO/CT2
      RETURN
      END
CCCCC
      SUBROUTINE VECP(Y,AMYS,YDOT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(12),YDOT(12)
      COMMON/MASS/EM(3),XM,EMU,DM(3),EPS2,EPS3
      COMMON/ENER/EVEC,HSTRE
      COMMON/FACT/SWT1,SWT2
      HSTR2=HSTRE*HSTRE
      ARO=Y(1)
      ARO2=ARO*ARO
      ARO3=ARO2*ARO
      PPHY=Y(9)
      ST=DSIN(Y(2))
      ST2=ST*ST
      ST3=ST2*ST
      CT=DCOS(Y(2))
      YDOT(1)=0.0D0
      YDOT(2)=0.0D0
      YDOT(3)=-SWT2*4.0D0*HSTRE/(AMYS*ARO2*ST2)
      YDOT(4)=0.0D0
      YDOT(5)=0.0D0
      YDOT(6)=0.0D0
      YDOT(7)=4.0D0*HSTR2*(SWT1-SWT2*2.0D0*PPHY/HSTRE)/
     1(AMYS*ARO3*ST2)
      YDOT(8)=4.0D0*HSTR2*CT*(SWT1-SWT2*2.0D0*PPHY/HSTRE)/
     1(AMYS*ARO2*ST3)
      YDOT(9)=0.0D0
      YDOT(10)=0.0D0
      YDOT(11)=0.0D0
      YDOT(12)=0.0D0
      EVEC=2.0D0*HSTR2*(SWT1-SWT2*2.0D0*PPHY/HSTRE)/
     1(AMYS*ARO2*ST2)
      RETURN
      END
CCCCC
      SUBROUTINE ROTVIB(Y,EROT,EVIB,RJA,AV,RLA)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/AIN/AA1,AA2,AA3
      COMMON/SPEC/BE,CE,AE,RE,AV1,D,BET
      COMMON/ENERGY/EKINA,EKINR  
      DIMENSION RJ(3),Y(12),RL(3)
      AM=AA1*AA2/(AA1+AA2) 
      AMY=AA3*(AA1+AA2)/(AA1+AA2+AA3)  
      RJ(1)=AM*(Y(8)*Y(12)-Y(9)*Y(11)) 
      RJ(2)=-AM*(Y(7)*Y(12)-Y(9)*Y(10))
      RJ(3)=AM*(Y(7)*Y(11)-Y(8)*Y(10)) 
      RL(1)=AMY*(Y(2)*Y(6)-Y(3)*Y(5))  
      RL(2)=-AMY*(Y(1)*Y(6)-Y(3)*Y(4)) 
      RL(3)=AMY*(Y(1)*Y(5)-Y(2)*Y(4))  
      RJA=DSQRT(RJ(1)**2+RJ(2)**2+RJ(3)**2)  
      RLA=DSQRT(RL(1)**2+RL(2)**2+RL(3)**2)  
      RA=DSQRT(Y(7)**2+Y(8)**2+Y(9)**2)
      ETOT=0.5*AM*(Y(10)**2+Y(11)**2+Y(12)**2)
      EKINA=ETOT
      EKINR=0.5*AMY*(Y(4)**2+Y(5)**2+Y(6)**2)
      VAB=D*(1.-DEXP(-BET*(RA-RE)))**2 
      ETOT=ETOT+VAB  
      RJA=15.74609416*RJA  
      RLA=15.74609416*RLA  
      RJA=RJA-0.5
      RLA=RLA-0.5
      EROTG=0.0
      AV=0.0
    2 BR=BE-AE*(AV+0.5)
      EROT=(RJA+0.5)**2*(BR-CE*(RJA+0.5)**2) 
      EVIB=ETOT-EROT 
      ARG=1.-EVIB/D  
C     IF (ARG.LT.0.0) WRITE(6,10) RJA,RLA,ETOT,EROT,EVIB 
   10 FORMAT(' ARG.LT.0.0 ',6E10.3)
      IF (ARG.LT.0.0) ARG=0.0
      AV=DSQRT(ARG)  
      AV=.5*(AV1*(1.-AV)-1.)
      IF (DABS(EROTG-EROT).LT.0.0001D0) GO TO 3
      EROTG=EROT
      GO TO 2  
    3 RETURN
      END
C TO DELF-COORDINATES IN D(12)
      SUBROUTINE CARDEL(Y,D)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/AMASS/AM1,AM2,AM3
c     common/rotat/omi(3,3)
      DIMENSION Y(12),W(6),AS(3),RS(3),A(3),R(3),D(12)
      DIMENSION PA(3),PR(3),OM(3,3),O(3,3)
      DIMENSION T(6,6),YG(12),TIN(6,6),Q(6)  
      DIMENSION OMI(3,3) 
      DATA PI/3.14159265359D0/
C CALC. OF TOTAL ANGULAR MOMENTUM
      AMR=AM3*AM2/(AM3+AM2)
      AM=AM1+AM2+AM3 
      AMY=AM1*(AM2+AM3)/AM 
      RJX=AMR*(Y(8)*Y(12)-Y(9)*Y(11))  
      RJY=AMR*(Y(9)*Y(10)-Y(7)*Y(12))  
      RJZ=AMR*(Y(7)*Y(11)-Y(8)*Y(10))  
      RLX=AMY*(Y(2)*Y(6)-Y(3)*Y(5))
      RLY=AMY*(Y(3)*Y(4)-Y(1)*Y(6))
      RLZ=AMY*(Y(1)*Y(5)-Y(2)*Y(4))
      RJ1=RJX+RLX
      RJ2=RJY+RLY
      RJ3=RJZ+RLZ
      RJT=DSQRT(RJ1*RJ1+RJ2*RJ2+RJ3*RJ3)
c     WRITE(6,15) RJ1,RJ2,RJ3,RJT
   15 FORMAT(' TOTAL MOMENTA ',4E10.3) 
      PHI=DATAN(RJ2/RJ1)
      IF (RJ1.LT.0.0) PHI=PHI+PI 
      TE=DACOS(RJ3/RJT)
      SIP=DSIN(PHI)  
      COP=DCOS(PHI)  
      SIT=DSIN(TE)
      COT=DCOS(TE)
      OM(1,1)=-SIP
      OM(1,2)=COP
      OM(1,3)=0.0
      OM(2,1)=-COT*COP
      OM(2,2)=-COT*SIP
      OM(2,3)=SIT
      OM(3,1)=SIT*COP
      OM(3,2)=SIT*SIP
      OM(3,3)=COT
      DO 10 I=1,3
      SU=0.0
      SU1=0.0  
      SU2=0.0  
      SU3=0.0  
      DO 11 J=1,3
      omi(j,i)=om(i,j)
      SU=SU+OM(I,J)*Y(J)
      SU1=SU1+OM(I,J)*Y(J+3)
      SU2=SU2+OM(I,J)*Y(J+6)
   11 SU3=SU3+OM(I,J)*Y(J+9)
      YG(I)=SU 
      YG(I+3)=SU1
      YG(I+6)=SU2
   10 YG(I+9)=SU3
      DO 12 I=1,12
   12 Y(I)=YG(I)
      AM=AM1+AM2+AM3 
      AMY=DSQRT(AM1*AM2*AM3/AM)  
      D1=DSQRT(AM1*(1.-AM1/AM)/AMY)
      A(1)=Y(7)/D1
      A(2)=Y(8)/D1
      A(3)=Y(9)/D1
      R(1)=Y(1)*D1
      R(2)=Y(2)*D1
      R(3)=Y(3)*D1
      PA(1)=Y(10)/D1 
      PA(2)=Y(11)/D1 
      PA(3)=Y(12)/D1 
      PR(1)=Y(4)*D1  
      PR(2)=Y(5)*D1  
      PR(3)=Y(6)*D1  
C A(1)=RX  R(1)=RRX  PA(1)=D RX/DT  PR(1)=D RRX/DT 
      AX=0.5*(A(2)*R(3)-A(3)*R(2))
      AY=0.5*(A(3)*R(1)-A(1)*R(3))
      AZ=0.5*(A(1)*R(2)-A(2)*R(1))
      AN=DSQRT(AX*AX+AY*AY+AZ*AZ)
      BETA=DACOS(AZ/AN)
      ALFA=DATAN(AY/AX)
      IF (AX.LT.0.0D0) ALFA=ALFA+PI
      SA=DSIN(ALFA)  
      CA=DCOS(ALFA)  
      CB=DCOS(BETA)  
      SB=DSIN(BETA)  
      OM(1,1)=-SA
      OM(1,2)=CA
      OM(1,3)=0.0
      OM(2,1)=-CB*CA 
      OM(2,2)=-CB*SA 
      OM(2,3)=SB
      OM(3,1)=SB*CA  
      OM(3,2)=SB*SA  
      OM(3,3)=CB
c      write(6,100) (a(i),r(i),i=1,3)
      DO 1 I=1,3
      SU=0.0
      SU1=0.0  
      DO 2 J=1,3
      SU=SU+OM(I,J)*A(J)
    2 SU1=SU1+OM(I,J)*R(J) 
      AS(I)=SU 
    1 RS(I)=SU1
c     WRITE(6,100) (AS(i),RS(i),i=1,3)
 100  FORMAT(1X,10E10.3)
      SS=AS(1)**2+RS(1)**2-AS(2)**2-RS(2)**2 
      TS=2.*(AS(1)*AS(2)+RS(1)*RS(2))  
      GA=0.5*DATAN(TS/SS)  
      IF (SS.LT.0.0D0) GA=GA+0.5*PI
      RHO=DSQRT(AS(1)**2+AS(2)**2+RS(1)**2+RS(2)**2)
      D(4)=ALFA
      D(5)=BETA
      D(6)=GA  
      CG=DCOS(GA)
      SG=DSIN(GA)
      SA=DSIN(ALFA)  
      CA=DCOS(ALFA)  
      CB=DCOS(BETA)  
      SB=DSIN(BETA)  
      O(1,1)=-CG*SA-CB*CA*SG
      O(1,2)=CG*CA-CB*SA*SG
      O(1,3)=SG*SB
      O(2,1)=SG*SA-CB*CA*CG
      O(2,2)=-SG*CA-CB*SA*CG
      O(2,3)=CG*SB
      O(3,1)=SB*CA
      O(3,2)=SB*SA
      O(3,3)=CB
      DO 7 I=1,3
      SU=0.0
      SU1=0.0  
c      write(6,100) (o(i,j),j=1,3)
      DO 8 J=1,3
      SU=SU+O(I,J)*A(J)
    8 SU1=SU1+O(I,J)*R(J)  
      AS(I)=SU 
    7 RS(I)=SU1
c      WRITE(6,302) (AS(I),I=1,3) 
c      WRITE(6,302) (RS(I),I=1,3) 
      TE=DATAN(RS(2)/AS(1))
      AS1=AS(1)
C     IF (AS1.LT.0.0D0) TE=TE+PI 
      PHI=DATAN(RS(1)/AS(1))
      IF (AS1.LT.0.0D0) PHI=PHI+PI
c      IF (PHI.GT.0.0) PHI=PHI-PI 
C     WRITE(6,101) TE,PHI  
  101 FORMAT(' TE,PHI CARDEL=',2E10.3) 
      CT=DCOS(TE)
      ST=DSIN(TE)
      CP=DCOS(PHI)
      SP=DSIN(PHI)
      RCC=RHO*CT*CP  
      RSC=RHO*ST*CP  
      RCS=RHO*CT*SP  
      RSS=RHO*ST*SP  
      AX=O(1,2)*A(1)-O(1,1)*A(2) 
      AY=O(2,2)*A(1)-O(2,1)*A(2) 
      AZ=O(3,2)*A(1)-O(3,1)*A(2) 
      AXM=O(1,2)*R(1)-O(1,1)*R(2)
      AYM=O(2,2)*R(1)-O(2,1)*R(2)
      AZM=O(3,2)*R(1)-O(3,1)*R(2)
      BZ=-CB*CA*A(1)-CB*SA*A(2)+SB*A(3)
      BZM=-CB*CA*R(1)-CB*SA*R(2)+SB*R(3)
      DO 4 I=1,6
      tin(i,i)=1.0d0
      DO 4 J=1,6
    4 T(I,J)=0.0D0
      T(1,1)=RCC
      T(1,2)=-RSC
      T(1,3)=-RCS
      T(1,4)=AX
      T(1,6)=RSS
      T(2,1)=-RSS
      T(2,2)=-RCS
      T(2,3)=-RSC
      T(2,4)=AY
      T(2,6)=RCC
      T(4,1)=RCS
      T(4,2)=-RSS
      T(4,3)=RCC
      T(4,4)=AXM
      T(4,6)=-RSC
      T(3,4)=AZ
      T(3,5)=BZ
      T(5,1)=RSC
      T(5,2)=RCC
      T(5,3)=-RSS
      T(5,4)=AYM
      T(5,6)=RCS
      T(6,4)=AZM
      T(6,5)=BZM
      TE=0.5*PI-2.*TE
C     PHI=2.*PHI
C     PHI=0.5*PI-2.*PHI
      PHI=-2.*PHI
      D(1)=RHO 
      D(2)=TE  
      D(3)=PHI 
      IFAIL=0  
      DO 22 I=1,3
      Q1=0.0
      Q2=0.0
      DO 3 K=1,3
      Q1=Q1+O(I,K)*PA(K)
    3 Q2=Q2+O(I,K)*PR(K)
      Q(I)=Q1  
   22 Q(I+3)=Q2
C     PRINT 300
  300 FORMAT(' T-MATRIX')  
C     DO 301 I=1,6
C 301 WRITE(6,302) (T(I,K),K=1,6),Q(I) 
  302 FORMAT(1X,7E10.3)
      M=6
c      CALL F01AAF(T,6,M,TIN,6,W,IFAIL) 
      call gaussj(t,6,m,tin,6,m)
      DO 5 I=1,6
      SU=0.0
      DO 6 K=1,6
      tin(i,k)=t(i,k)
    6 SU=SU+TIN(I,K)*Q(K)  
    5 D(I+6)=SU
      D(8)=-2.*D(8)  
C     D(9)=2.*D(9)
      D(9)=-2.*D(9)  
C     DO 200 I=1,3
C 200 WRITE(6,100) (O(I,J),J=1,3)
C     WRITE(6,100) (D(I),I=1,12) 
C TRANSFORMATION TO MOMENTA
      ARHO=AMY*D(1)*D(1)
      PAR=1.-DSIN(D(2))*DCOS(2.*D(6))  
      D(7)=AMY*D(1)*D(7)
      D(8)=ARHO*0.25*D(8)  
      CG2=DCOS(D(2)) 
      D(10)=0.5*D(10)*ARHO*CG2*CG2/PAR 
      D(12)=0.5*ARHO*CG2*CG2*(2.*D(12)+CG2*D(9))/(CG2*CG2-2.*PAR)
      D(9)=0.25*ARHO*D(9)*(1.-CG2*CG2)+0.5*D(12)*CG2
C COMMENT FOLLOWING 3 STATEMENTS IF PGAMMA NEQ JTOT
c      D(10)=RJT
c      D(11)=0.0 
c      D(12)=RJT*DCOS(D(5)) 
C     WRITE(6,303) D(9),D(10),D(12)
  303 FORMAT(' PP,PA,PG=',3E10.3)
      RETURN
      END
      SUBROUTINE D02AHF(X, Y, E, T, N, H0, H, AUX, F, YP, D, DP, N1, Q,
C     D02AHF INTEGRATES THE DIFF EQUATIONS DEFINED 
C     BY AUX OVER A RANGE, USING VARIABLE-ORDER ADAMS
C     MARK 4.5 REVISED
     *IFAIL)
      INTEGER T, N, N1, Q, IFAIL, P01AAF, Q0, ST, ER, ERR, SK, SMAX, U,
     *X02BBF, M, I, K, IW, JZ, J, IQ, J1, K1, L, I1
      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION X, H0, H, EPS, S, P, R, X0, WS, SS, PP, Y(N), E(N
     *), 
     * X02AAF, 
     * F(N), YP(N), D(N1,15), DP(N1,15), GMA(14), ETA(23), C(14,14)  
C  
C     SET COEFFICIENTS AND STEP-LENGTH 
      external aux
      DATA SRNAME /8H D02AHF /
C     SMALLEST REAL SUCH THAT 1+EPS>1  
c      EPS = X02AAF(EPS)
      eps=.5e-38
C     LARGEST INTEGER
c      SMAX = X02BBF(EPS)
      smax=2147483647
      M = 0
      IF (H0.EQ.0.0D0) GO TO 1180
      S = dMAX1(dABS(X),dABS(X+H0))
      IF (ABS(H0).LT.EPS*S) GO TO 1180 
      P = 1.0 D0
      IF (IABS(T).EQ.2) P = 0.0 D0
      R = 1.0 D0
      IF (IABS(T).EQ.3) R = 0.0 D0
      X0 = X + H0
      SK = 1
      GMA(1) = 1.0 D0
      GMA(2) = 0.5 D0
      GMA(3) = 5.0D0/12.0 D0
      GMA(4) = 0.375 D0
      GMA(5) = 251.0D0/720.0 D0  
      GMA(6) = 475.0D0/1440.0 D0 
      GMA(7) = 19087.0D0/60480.0 D0
      GMA(8) = 36799.0D0/120960.0 D0
      GMA(9) = 1070017.0D0/3628800.0 D0
      GMA(10) = 2082753.0D0/7257600.0 D0
      DO 40 I=10,13  
       S = 0.0 D0
       DO 20 K=1,I
       WS = K + 1
       IW = I - K + 1
       S = S + GMA(IW)/WS  
   20  CONTINUE
       GMA(I+1) = 1.0 D0- S
   40 CONTINUE 
      IF (H.NE.0.0 D0.AND. ABS(H/H0).LE.1.0D0) GO TO 80  
      IF (T.GE.0) GO TO 60 
      M = 2
      GO TO 1180
   60 H = H0
   80 S = ABS(H0/H) + 0.1 D0
      IF (S.LE.DFLOAT(SMAX/2)) GO TO 100
      M = 1
      GO TO 1180
  100 ST = INT(S)
      K = 1
      IF (T.GE.0) GO TO 120
      Q0 = Q/100
      Q = Q - 100*Q0 
      GO TO 920
  120 IF (K.GE.ST) GO TO 140
      K = 2*K  
      GO TO 120
  140 ST = K
      H = H0/DFLOAT(ST)
      T = -T
      Q = 1
      Q0 = -1  
      CALL AUX(F, Y, X)
      DO 160 I=1,N
       D(I,1) = F(I) 
       D(I,2) = 0.0 D0
C  
C     ADVANCE - TAKE STEP, TEST FOR ORDER AND STEP CHANGE
  160 CONTINUE 
  180 JZ = 1
      GO TO 1200
  200 DO 220 I=1,N
       WS = E(I)*(P*ABS(YP(I))+R)
       WS = 10.0D0*H*DP(I,Q+1)*(GMA(Q+1)-GMA(Q))/WS
       IF (ABS(WS).GT.1.0D0) GO TO 240 
  220 CONTINUE 
      JZ = 4
      GO TO 1300
  240 H = H/2.0 D0
      IF (ST.LE.SMAX/2) GO TO 260
      M = 1
      GO TO 1180
  260 ST = 2*ST
      IF (Q.LE.1) GO TO 180
      Q0 = 0
      DO 380 J=1,N
       IQ = Q - 1
       DO 280 K=1,IQ 
       J1 = J  
       C(1,K) = D(J1,K+1)*(0.25D0**K)  
  280  CONTINUE
       I = Q - 2
  300  IF (I.EQ.0) GO TO 340
       K = Q - 2
  320  C(1,K) = C(1,K) + C(1,K+1)
       K = K - 1
       IF (K.GE.I) GO TO 320
       I = I - 1
       GO TO 300
  340  DO 360 K=1,IQ 
       J1 = J  
       D(J1,K+1) = C(1,K)*(2.0D0**K)
  360  CONTINUE
  380 CONTINUE 
      IF (Q.LE.3) GO TO 180
      C(1,1) = 0.5 D0
      DO 440 J=1,Q
       C(J+1,J+1) = 0.5D0*C(J,J) 
       I = J
  400  IF (I.LT.2) GO TO 420
       C(I,J+1) = 0.5D0*(C(I+1,J+1)+C(I-1,J))
       I = I - 1
       GO TO 400
  420  C(1,J+1) = 0.5D0*C(2,J+1) 
  440 CONTINUE 
  460 JZ = 2
      GO TO 1200
  480 U = Q
      IF (Q.EQ.1) GO TO 580
      DO 560 J=1,N
       ERR = 0 
       IW = U - 1
       DO 540 K=1,IW 
       J1 = J  
       ER = 1  
       IF (ABS(DP(J1,K+2)).GT.ABS(DP(J1,K+1))) ER = 0
       IF (ER.EQ.1) GO TO 520
       IF (ERR.EQ.0) GO TO 500
       IF (Q.GT.K-1) Q = K - 1
       GO TO 560
  500  ERR = 1 
       GO TO 540
  520  IF (ERR.EQ.1) ERR = 0
  540  CONTINUE
  560 CONTINUE 
      IF (Q.EQ.U) GO TO 580
      IF (SK.EQ.0) SK = 1  
      GO TO 180
  580 SK = 1 - SK
      IF (Q.EQ.1) GO TO 840
      IF (SK.EQ.1) GO TO 680
      S = 0.0 D0
      DO 600 K=1,Q
       S = S + C(K,Q)
  600 CONTINUE 
      IQ = Q - 1
      DO 620 K=1,IQ  
       ETA(K) = C(K,Q)/S
  620 CONTINUE 
      DO 660 J=1,N
       DO 640 K=1,IQ 
       J1 = J  
       D(J1,K+1) = D(J1,K+1) + ETA(K)*DP(J1,Q+1)
  640  CONTINUE
  660 CONTINUE 
      GO TO 840
  680 SS = 0.0 D0
      IQ = Q + 1
      DO 700 J=1,IQ  
       SS = SS + C(J,Q+1)  
  700 CONTINUE 
      DO 740 K=2,Q
       PP = 0.0 D0
       S = 0.0 D0
       IW = K - 1
       DO 720 J=1,IW 
       PP = PP + ETA(J)
       S = S + C(J,Q+1)
  720  CONTINUE
       C(K,1) = PP*SS - S  
  740 CONTINUE 
      IF (Q.LE.2) GO TO 840
      S = 0.0 D0
      DO 760 K=2,Q
       S = S + C(K,1)
  760 CONTINUE 
      S = S + C(Q+1,Q+1)
      IQ = Q - 1
      DO 780 K=2,IQ  
       IW = Q - 2 + K
       ETA(IW) = C(K,1)/S  
  780 CONTINUE 
      DO 820 J=1,N
       DO 800 K=2,IQ 
       J1 = J  
       IW = Q - 2 + K
       D(J1,K+1) = D(J1,K+1) + ETA(IW)*(DP(J1,Q+1)-D(J1,Q+1))  
  800  CONTINUE
  820 CONTINUE 
  840 DO 880 K1=1,N  
       S = 0.0 D0
       DO 860 J=1,Q  
       K = K1  
       J1 = Q + 1 - J
       S = S + GMA(J1)*D(K,J1)
  860  CONTINUE
       YP(K) = Y(K) + H*S  
  880 CONTINUE 
      JZ = 3
      GO TO 1300
  900 IF (SK.EQ.1) GO TO 180
C  
C     ENTRY POINT AFTER FIRST STEP
      GO TO 460
  920 IF (Q0.LT.Q) GO TO 180
      L = Q - 1
      IF (L.EQ.0) L = 1
      IF (Q.EQ.13) U = 13  
      IF (Q.NE.13) U = Q + 1
      S = 0.0 D0
      DO 940 I=1,N
       S = S + ABS(D(I,L+1)*(GMA(L+1)-GMA(L)))
  940 CONTINUE 
      S = S*(0.5D0**(Q-L+1))
      J = L + 1
  960 IF (J.GT.U) GO TO 1020
      SS = 0.0 D0
      DO 980 I=1,N
       SS = SS + ABS(D(I,J+1)*(GMA(J+1)-GMA(J)))
  980 CONTINUE 
      SS = SS*(0.5D0**(Q-J+1))
      IF (S.LE.SS) GO TO 1000
      S = SS
      L = J
 1000 J = J + 1
      GO TO 960
 1020 IF (L.GT.Q) Q0 = 0
      Q = L
      DO 1040 I=1,N  
       WS = E(I)*(P*ABS(Y(I))+R) 
       WS = 10.0D0*H*D(I,Q+1)*(GMA(Q+1)-GMA(Q))/WS 
       IF (ABS(WS).GT.0.5D0**(Q+1)) GO TO 180
 1040 CONTINUE 
      IW = ST/2
      IF (ST.NE.2*IW) GO TO 180  
      H = 2.0D0*H
      ST = ST/2
      IF (Q.LE.1) GO TO 180
      Q0 = 0
      DO 1160 J=1,N  
       IQ = Q - 1
       DO 1060 K=1,IQ
       J1 = J  
       C(1,K) = D(J1,K+1)*(0.5D0**K)
 1060  CONTINUE
       I = 1
 1080  IF (I.GE.Q-1) GO TO 1120  
       IW = Q - 2
       DO 1100 K=I,IW
       C(1,K) = C(1,K) - C(1,K+1)
 1100  CONTINUE
       I = I + 1
       GO TO 1080
 1120  DO 1140 K=1,IQ
       J1 = J  
       D(J1,K+1) = C(1,K)*(4.0D0**K)
 1140  CONTINUE
 1160 CONTINUE 
      GO TO 180
 1180 IF (M.EQ.0) IFAIL = M
      IF (M.NE.0) IFAIL = IFAIL+M
      IF (IFAIL.EQ.10) WRITE(6,2222) M
      IF (IFAIL.EQ.10) STOP
 2222 FORMAT(' IFAIL=',I3)
C  
C     D02AHZ - PREDICTS NEW VALUE
      RETURN
 1200 WS = ST - 1
      X = X0 - WS*H  
      DO 1240 I=1,N  
       S = 0.0 D0
       DO 1220 J=1,Q 
       I1 = I  
       J1 = Q + 1 - J
       S = S + GMA(J1)*D(I1,J1)  
 1220  CONTINUE
       YP(I) = Y(I) + H*S  
 1240 CONTINUE 
      CALL AUX(F, YP, X)
      DO 1280 I=1,N  
       DP(I,1) = F(I)
       DO 1260 J=1,Q 
       I1 = I  
       DP(I1,J+1) = DP(I1,J) - D(I1,J) 
 1260  CONTINUE
 1280 CONTINUE 
      IF (JZ.EQ.1) GO TO 200
C  
C     D02AHY - CORRECTOR ROUTINE 
      GO TO 480
 1300 DO 1320 I=1,N  
       Y(I) = YP(I) + H*GMA(Q+1)*DP(I,Q+1)
 1320 CONTINUE 
      CALL AUX(F, Y, X)
      DO 1380 I=1,N  
       DP(I,1) = F(I)
       IQ = Q + 1
       DO 1340 J=1,IQ
       I1 = I  
       DP(I1,J+1) = DP(I1,J) - D(I1,J) 
 1340  CONTINUE
       IQ = Q + 2
       DO 1360 J=1,IQ
       I1 = I  
       D(I1,J) = DP(I1,J)  
 1360  CONTINUE
 1380 CONTINUE 
      IF (Q0.LT.14) Q0 = Q0 + 1  
      ST = ST - 1
      IF (ST.NE.0) GO TO 1400
      Q = Q + 100*Q0 
      GO TO 1180
 1400 IF (JZ.EQ.3) GO TO 900
      GO TO 920
      END
      DOUBLE PRECISION FUNCTION RANDOM(L)
      IMPLICIT REAL*8 (A-H,O-Z)  
      COMMON /ARAND/IX,IY,IZ
C  EERSTE GENERATIE  
      IX=171*MOD(IX,177)-2*(IX/177)
      IF(IX.LT.0) IX=IX+30269
C  TWEEDE GENERATIE  
      IY=172*MOD(IY,176)-35*(IY/176)
      IF(IY.LT.0) IY=IY+30307
C  DERDE GENERATIE
      IZ=170*MOD(IZ,178)-63*(IZ/178)
      IF(IZ.LT.0) IZ=IZ+30323
C  COMBINATIE  
      RANDOM=AMOD(FLOAT(IX)/30269.0
     *+FLOAT(IY)/30307.0+FLOAT(IZ)/30323.0,1.0)
      RETURN
      END
c documentation numerical recipies p. 28
c a is an input matrix n times n stored in array with
c physical dimensions np. b is input matrix n by m stored in
c array with physical dimensions np by mp. On output a is the
c inverse matrix and b is replaced by solution vectors
      subroutine gaussj(a,n,np,b,m,mp)
      implicit real*8(a-h,o-z)
      parameter (nmax=50)
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
      do 11 j=1,n
      ipiv(j)=0
   11 continue
      do 22 i=1,n
      big=0.0
      do 13 j=1,n
      if (ipiv(j).ne.1)then
      do 12 k=1,n
      if (ipiv(k).eq.0) then
      if (dabs(a(j,k)).ge.big) then
      big=dabs(a(j,k))
      irow=j
      icol=k
      endif
      else if (ipiv(k).gt.1) then
      pause 'Singular matrix'
      endif
   12 continue
      endif
   13 continue
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
      do 14 l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
   14 continue
      do 15 l=1,m
        dum=b(irow,l)
        b(irow,l)=b(icol,l)
        b(icol,l)=dum
   15 continue
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.d0) pause 'Singular matrix'
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.d0
      do 16 l=1,n
      a(icol,l)=a(icol,l)*pivinv
   16 continue
      do 17 l=1,m
      b(icol,l)=b(icol,l)*pivinv
   17 continue
      do 21 ll=1,n
      if (ll.ne.icol) then
      dum=a(ll,icol)
      a(ll,icol)=0.d0
      do 18 l=1,n
      a(ll,l)=a(ll,l)-a(icol,l)*dum
   18 continue
      do 19 l=1,m
      b(ll,l)=b(ll,l)-b(icol,l)*dum
   19 continue
      endif
   21 continue
   22 continue
      do 24 l=n,1,-1
        if (indxr(l).ne.indxc(l))then
        do 23 k=1,n
        dum=a(k,indxr(l))
        a(k,indxr(l))=a(k,indxc(l))
        a(k,indxc(l))=dum
   23 continue
      endif
   24 continue
      return
      end
      SUBROUTINE DELCAR(D,Y)
      IMPLICIT REAL*8(A-H,O-Z)
C TRANSFORMS FROM DELFS TO CARTESIAN COORDINATES
      COMMON/AMASS/AM1,AM2,AM3
      COMMON/MASS/EM(3),XM,EMU,DM(3),EPS2,EPS3
      COMMON/REA/ICH 
c     common/rotat/omi(3,3)
      DIMENSION Y(12),D(12)
      DIMENSION O(3,3),OA(3,3),OB(3,3),OG(3,3),A(3),B(3),C(3)  
      DIMENSION YP(12)
      DIMENSION E(3),OMA(3),OMR(3),R(3),F(3) 
      DIMENSION T(6,6),W(6),TIN(6,6),Q(6)
      DIMENSION OM(3,3),OMI(3,3)
      DATA PI/3.14159265359D0/
C     WRITE(6,200) (D(I),I=1,12) 
  200 FORMAT('D-VECTOR=',6F10.6)
      SIP=DSIN(D(3))
      COP=DCOS(D(3))
      SIT=DSIN(D(2))
      COT=DCOS(D(2))
      OM(1,1)=-SIP
      OM(1,2)=COP
      OM(1,3)=0.0
      OM(2,1)=-COT*COP
      OM(2,2)=-COT*SIP
      OM(2,3)=SIT
      OM(3,1)=SIT*COP
      OM(3,2)=SIT*SIP
      OM(3,3)=COT
      DO 57 I=1,3
      DO 58 J=1,3
      OMI(J,I)=OM(I,J)
  58  CONTINUE
  57  CONTINUE  
      AM=AM1+AM2+AM3 
      AMY=DSQRT(AM1*AM2*AM3/AM)  
      D1=DM(1) 
      ARHO=AMY*D(1)**2
      AZ=ARHO  
      ST=DSIN(D(2))  
      STT=ST
      ST2=ST*ST
      CT=DCOS(D(2))  
      CTT=CT
      CT2=CT*CT
      PP=D(9)  
      CG=DCOS(D(6))  
      CG2=DCOS(2.*D(6))
      SG=DSIN(D(6))  
      CA=DCOS(D(4))  
      SA=DSIN(D(4))  
      CB=DCOS(D(5))  
      SB=DSIN(D(5))  
      O(1,1)=-CG*SA-CB*CA*SG
      O(1,2)=SG*SA-CB*CA*CG
      O(1,3)=SB*CA
      O(2,1)=CG*CA-CB*SA*SG
      O(2,2)=-SG*CA-CB*SA*CG
      O(2,3)=SB*SA
      O(3,1)=SG*SB
      O(3,2)=CG*SB
      O(3,3)=CB
      BETA=-DACOS(D(12)/D(10))
      AJX=-D(10)*SG*DSIN(BETA)
      AJY=-D(10)*CG*DSIN(BETA)
      AJZ=-D(10)*DCOS(BETA)
C     AJZ=-D(12)
      AX=0.5*ARHO*(1.-ST)  
      AY=0.5*ARHO*(1.+ST)  
      DO 42 I=1,6
      tin(i,i)=1.0d0
      DO 42 J=1,6
  42  T(I,J)=0.0
      T(1,1)=1.0D0
      T(2,2)=.25*D(1)*D(1) 
      T(3,3)=T(2,2)  
      T(4,4)=AX/AMY  
      T(5,5)=AY/AMY  
      T(6,6)=AZ/AMY  
      T(3,6)=0.5*D(1)*D(1)*CT
      T(6,3)=T(3,6)  
      M=6
      IFAIL=0  
c      CALL F01AAF(T,6,M,TIN,6,W,IFAIL) 
      call gaussj(t,6,m,tin,6,m)
      Q(1)=D(7)
      Q(2)=D(8)
      Q(3)=D(9)
      Q(4)=AJX 
      Q(5)=AJY 
      Q(6)=-AJZ
      DO 43 I=1,6
      SU=0.0
      DO 44 J=1,6
      tin(i,j)=t(i,j)
  44  SU=SU+TIN(I,J)*Q(J)  
      W(I)=SU  
   43 CONTINUE 
      TE=0.25*PI-0.5*D(2)  
      PHI=-0.5*D(3)  
C JOHNSON CHOICE FOR PHI
      CT=DCOS(TE)
      CP=DCOS(PHI)
      ST=DSIN(TE)
      SP=DSIN(PHI)
      A(1)=D(1)*CT*CP
      A(2)=-D(1)*ST*SP
      A(3)=0.0 
      R(1)=D(1)*CT*SP
      R(2)=D(1)*ST*CP
      R(3)=0.0 
c      write(6,500) (a(i),r(i),i=1,3)
  500 format(1x,10e10.3)
      DO 1 I=1,3
      SU=0.0
      SU1=0.0  
c      write(6,500) (o(i,k),k=1,3)
      DO 2 K=1,3
      SU=SU+O(I,K)*A(K)
    2 SU1=SU1+O(I,K)*R(K)  
      Y(I)=SU1/D1
    1 Y(I+6)=SU*D1
c      write(6,500) (y(i),y(i+6),i=1,3)
      W(1)=W(1)/D(1) 
      W(2)=-0.5*W(2) 
      W(3)=-0.5*W(3) 
      Y(10)=(A(1)*W(1)-R(2)*W(2)-R(1)*W(3)-A(2)*W(6))*D1/AMY
      Y(11)=(A(2)*W(1)-R(1)*W(2)-R(2)*W(3)+A(1)*W(6))*D1/AMY
      Y(12)=(A(2)*W(4)-A(1)*W(5))*D1/AMY
      Y(4)=(R(1)*W(1)+A(2)*W(2)+A(1)*W(3)-R(2)*W(6))/D1/AMY
      Y(5)=(R(2)*W(1)+A(1)*W(2)+A(2)*W(3)+R(1)*W(6))/D1/AMY
      Y(6)=(R(2)*W(4)-R(1)*W(5))/D1/AMY
c rotation of momenta
      do 100 k=1,3
      a(k)=y(k+3)
      r(k)=y(k+9)
  100 continue
      do 101 i=1,3
      su=0.0
      su1=0.0
      do 102 k=1,3
      su=su+o(i,k)*a(k)
  102 su1=su1+o(i,k)*r(k)
      y(i+3)=su
      y(i+9)=su1
  101 continue
      do 302 i=1,3
      su=0.0
      su1=0.0
      su2=0.0
      su3=0.0
      do 303 k=1,3
      su=su+omi(i,k)*y(k)
      su1=su1+omi(i,k)*y(k+3)
      su2=su2+omi(i,k)*y(k+6)
      su3=su3+omi(i,k)*y(k+9)
  303 continue
      e(i)=su
      f(i)=su1
      a(i)=su2
      r(i)=su3
  302 continue
      do 304 i=1,3
      y(i)=e(i)
      y(i+3)=f(i)
      y(i+6)=a(i)
      y(i+9)=r(i)
  304 continue
C     WRITE(6,100) TE,PHI  
c 100  FORMAT(' TE,PHI=',2E10.3)  
      EKR=0.5*D(7)*D(7)/AMY
      EIN=2.*(D(8)*D(8)+D(9)*D(9)/ST2)/ARHO  
      EGA=0.5*D(12)*D(12)/ARHO/ST2
      EGB=(D(10)**2-D(12)**2)/ARHO/CT2 
C SIGN DEPENDENT TERMS
      AA1=-2.*D(12)*D(9)*CTT/ARHO/ST2  
      AA2=-(D(10)**2-D(12)**2)*STT*DCOS(2.*D(6))/ARHO/CT2
      EKT=EKR+EIN+EGA+EGB+AA1+AA2
C     WRITE(6,104) EKR,EIN,EGA,EGB,EKT 
  104 FORMAT(' EKR ...=',5F10.6) 
C     WRITE(6,103) AA1,AA2 
  103 FORMAT(' AA1,AA2=',2F10.6) 
c  101 FORMAT(1X,7E10.3)
C TRANSFORMATION TO OTHER CHANNELS
      EP2=AM2/(AM2+AM3)
      EP3=1.-EP2
      IF (ICH.EQ.1) RETURN 
      IF (ICH.EQ.3) GO TO 33
      YP(1)=-EP3*Y(7)-(AM1*Y(1)+AM3*EP2*Y(7))/(AM1+AM3)  
      YP(2)=-EP3*Y(8)-(AM1*Y(2)+AM3*EP2*Y(8))/(AM1+AM3)  
      YP(3)=-EP3*Y(9)-(AM1*Y(3)+AM3*EP2*Y(9))/(AM1+AM3)  
      YP(4)=-EP3*Y(10)-(AM1*Y(4)+AM3*EP2*Y(10))/(AM1+AM3)
      YP(5)=-EP3*Y(11)-(AM1*Y(5)+AM3*EP2*Y(11))/(AM1+AM3)
      YP(6)=-EP3*Y(12)-(AM1*Y(6)+AM3*EP2*Y(12))/(AM1+AM3)
      YP(7)=-EP2*Y(7)+Y(1) 
      YP(8)=-EP2*Y(8)+Y(2) 
      YP(9)=-EP2*Y(9)+Y(3) 
      YP(10)=-EP2*Y(10)+Y(4)
      YP(11)=-EP2*Y(11)+Y(5)
      YP(12)=-EP2*Y(12)+Y(6)
      GO TO 35 
  33  YP(1)=EP2*Y(7) -(AM1*Y(1) -AM2*EP3*Y(7) )/(AM1+AM2)
      YP(2)=EP2*Y(8) -(AM1*Y(2) -AM2*EP3*Y(8) )/(AM1+AM2)
      YP(3)=EP2*Y(9) -(AM1*Y(3) -AM2*EP3*Y(9) )/(AM1+AM2)
      YP(4)=EP2*Y(10)-(AM1*Y(4) -AM2*EP3*Y(10))/(AM1+AM2)
      YP(5)=EP2*Y(11)-(AM1*Y(5) -AM2*EP3*Y(11))/(AM1+AM2)
      YP(6)=EP2*Y(12)-(AM1*Y(6) -AM2*EP3*Y(12))/(AM1+AM2)
      YP(7)=-EP3*Y(7)-Y(1) 
      YP(8)=-EP3*Y(8)-Y(2) 
      YP(9)=-EP3*Y(9)-Y(3) 
      YP(10)=-EP3*Y(10)-Y(4)
      YP(11)=-EP3*Y(11)-Y(5)
      YP(12)=-EP3*Y(12)-Y(6)
   35 CONTINUE 
      DO 34 I=1,12
  34  Y(I)=YP(I)
      RETURN
      END
CCCC
      SUBROUTINE DERIV_POT(R,DV)
      IMPLICIT REAL*8(A-H,O-Z)      
      PARAMETER (ND=4)
      DIMENSION R(3),RA(3,ND),RRA(3)
      DIMENSION ENP(ND),ENN(ND),DV(3)
      COMMON/COEFF/COFTWO(27),COFTHR(169)
      DX1=0.004D0
      DX2=0.004D0
      DX3=0.004D0
      AFAC1=672.0D0
      AFAC2=168.0D0
      AFAC3=32.0D0
      AFAC4=3.0D0
      AFAC=840.0D0
*
      DO I=1,ND
      RA(1,I)=R(1)+DFLOAT(I)*DX1
      RA(2,I)=R(2)
      RA(3,I)=R(3)
      RRA(1)=RA(1,I)
      RRA(2)=RA(2,I)
      RRA(3)=RA(3,I)
      CALL POTENTIAL(RRA,EN)
C     CALL POTH3(RRA,EN)
      ENP(I)=EN
      ENDDO
      DO I=1,ND
      RA(1,I)=R(1)-DFLOAT(I)*DX1
      RA(2,I)=R(2)
      RA(3,I)=R(3)
      RRA(1)=RA(1,I)
      RRA(2)=RA(2,I)
      RRA(3)=RA(3,I)
      CALL POTENTIAL(RRA,EN)
C     CALL POTH3(RRA,EN)
      ENN(I)=EN
      ENDDO
      DV(1)=(AFAC4*ENN(4)-AFAC3*ENN(3)+AFAC2*ENN(2)-AFAC1*ENN(1)
     1     +AFAC1*ENP(1)-AFAC2*ENP(2)+AFAC3*ENP(3)-AFAC4*ENP(4))
     2     /AFAC/DX1
*
      DO I=1,ND
      RA(1,I)=R(1)
      RA(2,I)=R(2)+DFLOAT(I)*DX2
      RA(3,I)=R(3)
      RRA(1)=RA(1,I)
      RRA(2)=RA(2,I)
      RRA(3)=RA(3,I)
      CALL POTENTIAL(RRA,EN)
C     CALL POTH3(RRA,EN)
      ENP(I)=EN
      ENDDO
      DO I=1,ND
      RA(1,I)=R(1)
      RA(2,I)=R(2)-DFLOAT(I)*DX2
      RA(3,I)=R(3)
      RRA(1)=RA(1,I)
      RRA(2)=RA(2,I)
      RRA(3)=RA(3,I)
      CALL POTENTIAL(RRA,EN)
C     CALL POTH3(RRA,EN)
      ENN(I)=EN
      ENDDO
      DV(2)=(AFAC4*ENN(4)-AFAC3*ENN(3)+AFAC2*ENN(2)-AFAC1*ENN(1)
     1     +AFAC1*ENP(1)-AFAC2*ENP(2)+AFAC3*ENP(3)-AFAC4*ENP(4))
     2     /AFAC/DX1
*
      DO I=1,ND
      RA(1,I)=R(1)
      RA(2,I)=R(2)
      RA(3,I)=R(3)+DFLOAT(I)*DX3
      RRA(1)=RA(1,I)
      RRA(2)=RA(2,I)
      RRA(3)=RA(3,I)
      CALL POTENTIAL(RRA,EN)
C     CALL POTH3(RRA,EN)
      ENP(I)=EN
      ENDDO
      DO I=1,ND
      RA(1,I)=R(1)
      RA(2,I)=R(2)
      RA(3,I)=R(3)-DFLOAT(I)*DX3
      RRA(1)=RA(1,I)
      RRA(2)=RA(2,I)
      RRA(3)=RA(3,I)
      CALL POTENTIAL(RRA,EN)
C     CALL POTH3(RRA,EN)
      ENN(I)=EN
      ENDDO
      DV(3)=(AFAC4*ENN(4)-AFAC3*ENN(3)+AFAC2*ENN(2)-AFAC1*ENN(1)
     1     +AFAC1*ENP(1)-AFAC2*ENP(2)+AFAC3*ENP(3)-AFAC4*ENP(4))
     2     /AFAC/DX1
*
      RETURN
      END
CCCC
      SUBROUTINE POTENTIAL(RR,EN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MATWO=27)
      PARAMETER (MATHR=169)
      PARAMETER (NC=3)
      DIMENSION RR(NC)
      COMMON/COEFF/COFTWO(MATWO),COFTHR(MATHR)
C     -------------
      PI=4.0D0*DATAN(1.0D0)
      AUANG=0.529170D0
C     -------------
      CALL TWBODY(RR,VPOT2)
*
      RR1=RR(1)
      RR2=RR(2)
      RR3=RR(3)
      CALL THREHF(RR1,RR2,RR3,V3EHF)
*
      EN=VPOT2+V3EHF
      EN=27.2107*EN
      EN=EN+26.550D0
C     WRITE(80,80)VPOT2,V3EHF,EN
C 80  FORMAT(1X,4F12.6)
      RETURN
      END
CCC
      SUBROUTINE TWBODY(RR,VPOT2)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NC=3)
      PARAMETER (NCOF=10)
      PARAMETER (MA=27)
      DIMENSION RR(NC),RR_REF2(NC)
      DIMENSION TCOF(0:NCOF-1),DCCOF(4:11)
      DIMENSION AN(4:11),BN(4:11)
      DIMENSION VEHF(NC),VDC(NC)
      COMMON/COEFF/ACOF(MA),COFTHR(169)
      ONBODY=-1.0D0
CCCCC
      DO I=0,NCOF-1
      TCOF(I)=ACOF(I+1)
      ENDDO
      II=4
      DO I=NCOF+1,NCOF+8
      DCCOF(II)=ACOF(I)
      II=II+1
      ENDDO
      ALPH0=ACOF(NCOF+9)
      ALPH1=ACOF(NCOF+10)
      BET0=ACOF(NCOF+11)
      BET1=ACOF(NCOF+12)
      RREF2=ACOF(NCOF+13)
      PRHO=ACOF(NCOF+14)
      GAMA0=ACOF(MA-2)
      GAMA1=ACOF(MA-1)
      GAMA2=ACOF(MA)
C     ------------
C     EHF TERM
C     ------------
      DO KP=1,NC
      RR_REF2(KP)=RR(KP)-RREF2
      ENDDO
      DO KP=1,NC
      GAMA=GAMA0*(1.0D0+GAMA1*DTANH(GAMA2*RR_REF2(KP)))
      EXPFUNC=EXP(-GAMA*RR_REF2(KP))
      FUNC=0.0D0
      DO KK=0,NCOF-1
      FUNC=FUNC+TCOF(KK)*RR_REF2(KP)**KK
      ENDDO
      VEHF(KP)=FUNC*EXPFUNC/RR(KP)
      ENDDO
C     ------------
C     DC TERM
C     ------------
      DO II=4,11
        AN(II)=ALPH0/(FLOAT(II))**ALPH1
        BN(II)=BET0*EXP(-BET1*FLOAT(II))
      ENDDO
*
      VPOT2=0.0D0
      DO KP=1,NC
      VDC(KP)=0.0D0
      DO JJ=4,11
      CHI=(1-EXP(-AN(JJ)*RR(KP)/PRHO
     1   -BN(JJ)*RR(KP)*RR(KP)/PRHO**2))**JJ
      VDC(KP)=VDC(KP)-CHI*DCCOF(JJ)/RR(KP)**JJ
      ENDDO
      VPOT2=VPOT2+VEHF(KP)+VDC(KP)
      ENDDO
      VPOT2=VPOT2+ONBODY
      RETURN
      END
CCC
      SUBROUTINE THREHF(RR1,RR2,RR3,V3EHF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (MATHR=169)
      DIMENSION FNC(MATHR)
      DIMENSION COF1(102),COF2(102)
      COMMON/COEFF/COFTWO(27),COFTHR(MATHR)
CCCCC
      B_MOR=0.5D0
      RREF3=1.40D0
*
      RTIL1=(1-EXP(-B_MOR*(RR1/RREF3-1)))/B_MOR
      RTIL2=(1-EXP(-B_MOR*(RR2/RREF3-1)))/B_MOR
      RTIL3=(1-EXP(-B_MOR*(RR3/RREF3-1)))/B_MOR
*
      QQ1=SQRT(1.0D0/3.0D0)*RTIL1+SQRT(1.0D0/3.0D0)*RTIL2
     $   +SQRT(1.0D0/3.0D0)*RTIL3
      QQ2=SQRT(1.0D0/2.0D0)*RTIL2-SQRT(1.0D0/2.0D0)*RTIL3
      QQ3=SQRT(2.0D0/3.0D0)*RTIL1-SQRT(1.0D0/6.0D0)*RTIL2
     $   -SQRT(1.0D0/6.0D0)*RTIL3
*
      GM1=QQ1
      GM2=QQ2**2+QQ3**2
      GM3=QQ3*(QQ3**2-3.0D0*QQ2**2)
*
      TT2=1.0D0/(1+EXP(0.5D0*(GM1-12.0D0)))
*
      IJK=1
      III=1
      IV=1
      DO IQ=0,12
      DO JQ=0,12,2
      DO KQ=0,12,3
*
      IF ((IQ+JQ+KQ).LE.10) THEN
      COF1(III)=COFTHR(IJK)
      FNC(IJK)=GM1**IQ*GM2**(JQ/2)*GM3**(KQ/3)*TT2
      IJK=IJK+1
      COF2(IV)=COFTHR(IJK)
      FNC(IJK)=SQRT(GM2)*GM1**IQ*GM2**(JQ/2)*GM3**(KQ/3)*TT2
CCCC
      III=III+1
      IV=IV+1
      IJK=IJK+1
      ENDIF
*
      IF ((IQ+JQ+KQ).GT.10.AND.(IQ+JQ+KQ).LE.12) THEN
      COF1(III)=COFTHR(IJK)
      COF2(IV)=0.0D0
      FNC(IJK)=GM1**IQ*GM2**(JQ/2)*GM3**(KQ/3)*TT2
CCCC
      III=III+1
      IV=IV+1
      IJK=IJK+1
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      V3EHF=0.0D0
*
      DO IC=1,MATHR
      V3EHF=V3EHF+COFTHR(IC)*FNC(IC)
      ENDDO
      RETURN
      END
CCC
      SUBROUTINE POTH3(R,EN)  
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(3),DV(3),RR(3) 
      DIMENSION XX(8)
      DIMENSION DVA(3),DA(3),DQ(3),DVT(3)
      DIMENSION VT(3),VS(3),A(3),Q(3),DVS(3) 
      DATA BETA1,BETA2,BETA3/0.52D0,0.052D0,0.79D0/
      DATA A2,A3,A4,A5/.0012646477D0,-.0001585792D0,.0000079707D0,
     1-.0000001151D0/
      DATA B11,B12,B22,B23/3.0231771503D0,-1.08935219D0,1.7732141742D0,
     1-2.0979468223D0/
      DATA B24,B31,B32,B41,B42/-3.978850217D0,.4908116374D0,
     1-.8718696387D0,.1612118092D0,-.1273731045D0/ 
      DATA B51,B52/-13.3599568553D0,.9877930913D0/ 
      DATA ALF5/0.0035D0/  
      DATA C6,C8/6.89992032D0,219.9997304D0/ 
      DATA XX/2.80465D0,2.86869D0,1.93426D0,-0.85360D0,0.24768D0
     1,-0.02013D0,2.80521D0,0.01927D0/ 
      DO 1 I=1,3
    1 RR(I) = R(I) - 1.40059D0
      DO 2 I=1,3
C     IF (R(I) .LE. 7.0D0) THEN  
C     VS(I) = -0.174475D0*(1.0+XX(1)*RR(I)+XX(2)*RR(I)**2+
C    +    XX(3)*RR(I)**3-XX(4)*RR(I)**4+XX(5)*RR(I)**5+XX(6)*RR(I)**6
C    1+XX(8)*RR(I)**7)*
C    +        DEXP(-XX(7)*RR(I)) 
      EX7=DEXP(-XX(7)*RR(I))
      VS(I)=-0.174475D0*EX7*(RR(I)*(RR(I)*(RR(I)*(RR(I)*(RR(I) 
     1*(RR(I)*(RR(I)*XX(8)+XX(6))+XX(5))-XX(4))+XX(3))+XX(2))+ 
     1XX(1))+1.0D0)  
C     DVS(I) = -0.174475D0*(XX(1)+2.0D0*XX(2)*RR(I)+3.0D0*XX(3)*
C    +         RR(I)**2-4.0D0*XX(4)*RR(I)**3+5.0D0*XX(5)*RR(I)**4
C    1+6.0D0*XX(6)*RR(I)**5+7.0D0*XX(8)*RR(I)**6)* 
C    +         DEXP(-XX(7)*RR(I)) - XX(7)*VS(I)
      DVS(I)=-.174475D0*EX7*(RR(I)*(RR(I)*(RR(I)*(RR(I)*(RR(I) 
     1*(RR(I)*7.D0*XX(8)+6.D0*XX(6))+5.D0*XX(5))-4.D0*XX(4))
     1+3.D0*XX(3))+2.D0*XX(2))+XX(1)) -XX(7)*VS(I) 
C     ELSE
C     VS(I) = -6.49903*R(I)**(-6)-124.3991*R(I)**(-8)-
C    +        3285.828*R(I)**(-10)-0.8205*R(I)**2*DSQRT(R(I))* 
C    +        DEXP(-2.0*R(I))
C     DVS(I) = 6.0D0*6.49903*R(I)**(-7)+8.0D0*124.3991*R(I)**(-9)+
C    +         32858.28*R(I)**(-11)-5.0D0*0.8205*R(I)*DSQRT(R(I))/
C    +         2.0D0*DEXP(-2.0D0*R(I))+2.0D0*0.8205*R(I)**2*DSQRT(R(I))*
C    +         DEXP(-2.0D0*R(I)) 
C     ENDIF
      FSW=DEXP(-0.011*(R(I)-10.0D0)**4)
      R2=R(I)*R(I)
      R6=R2*R2*R2
      R8=R6*R2 
      IF (R(I).GT.10.0) FSW=1.0D0
      DFSW=-0.044*(R(I)-10.0D0)**3*FSW 
      IF (R(I).GT.10.0) DFSW=0.0D0
      VLR=-C6/R6-C8/R8
      VS(I)=VS(I)+FSW*VLR  
      R7=R6*R(I)
      R9=R7*R2 
      DVS(I)=DVS(I)+DFSW*VLR+FSW*(6.D0*C6/R7+8.D0*C8/R9) 
   2  CONTINUE 
      DO 3 I=1,3
      VT(I)=-1.2148730613*(-1.514663474+R(I)-1.46*R(I)**2)*
     1DEXP(-2.088442*R(I)) 
      DVT(I)=-2.088442*VT(I)-1.2148730613*(1.D0-2.D0*1.46D0*R(I))*
     1DEXP(-2.088442*R(I)) 
    3 CONTINUE 
      DO 4 I=1,3
      Q(I) = (VS(I) + VT(I))/2.0D0
      A(I) = (VS(I) - VT(I))/2.0D0
    4 CONTINUE 
      D = 0.5D0*((A(1)-A(2))**2+(A(2)-A(3))**2+(A(3)-A(1))**2) 
      VLEPS = Q(1) + Q(2) + Q(3) - DSQRT(D)  
      F1 = 1.0D0 + (R(1)**2-R(2)**2-R(3)**2)/(2.0*R(2)*R(3))+  
     +           (R(2)**2-R(3)**2-R(1)**2)/(2.0*R(3)*R(1))+
     +           (R(3)**2-R(1)**2-R(2)**2)/(2.0*R(1)*R(2))
      A1=DABS((R(1)-R(2))*(R(2)-R(3))*(R(3)-R(1))) 
C     VA=0.0012646477*A1**2-0.0001585792*A1**3+0.0000079707*A1**4
C    1-0.0000001151*A1**5  
      VA=A1*A1*(A1*(A1*(A1*A5+A4)+A3)+A2)
      RS=R(1)+R(2)+R(3)
      VA=VA*DEXP(-ALF5*RS**3)
      F2=1./R(1)+1./R(2)+1./R(3) 
      F3=(R(2)-R(1))**2+(R(3)-R(1))**2+(R(3)-R(2))**2
C     VB1=F1*(3.0231771503-1.08935219*RS)*DEXP(-0.52*RS) 
      EXS1=DEXP(-BETA1*RS) 
      VB1=F1*(B11+B12*RS)*EXS1
C     VB2=1.7732141742*F1**2-2.0979468223*F1**3-3.978850217*F1**4
      VB2=F1*F1*(F1*(F1*B24+B23)+B22)  
      EXS2=DEXP(-BETA2*RS*RS)
      VB2=VB2*EXS2
C     VB2=VB2*DEXP(-0.052*RS*RS) 
C     VB3=F2*(0.4908116374*F1*DEXP(-0.52*RS)-0.8718696387*F1**2*
C    1DEXP(-0.052*RS*RS))  
      VB3=F2*F1*(B31*EXS1+B32*EXS2*F1) 
C     VB4=F3*F1*(.1612118092*DEXP(-0.52*RS)-.1273731045*DEXP(-0.052* 
C    1RS*RS))  
      VB4=F3*F1*(B41*EXS1+B42*EXS2)
C     VB5=F1*(-13.3599568553+.9877930913*RS*RS)*DEXP(-.79*RS)  
      EXS3=DEXP(-BETA3*RS) 
      VB5=F1*(B51+B52*RS*RS)*EXS3
      EN=27.2107*(VLEPS+VA+VB1+VB2+VB3+VB4+VB5)
      DO 5 I=1,3
      DQ(I)=(DVS(I)+DVT(I))*0.5D0
    5 DA(I)=(DVS(I)-DVT(I))*0.5D0
C DERIVATIVE OF NON-LEPS TERMS
      B2=1.D0/R(1)+1.D0/R(2)+1.D0/R(3) 
      S23=R(1)*R(1)-R(2)*R(2)-R(3)*R(3)
      S13=R(2)*R(2)-R(1)*R(1)-R(3)*R(3)
      S12=R(3)*R(3)-R(1)*R(1)-R(2)*R(2)
      A23=R(2)*R(3)  
      A13=R(1)*R(3)  
      A12=R(1)*R(2)  
      B1=2.D0+S23/A23+S13/A13+S12/A12  
      B1=0.5*B1
      DADR1=0.0D0
      DADR2=0.0D0
      DADR3=0.0D0
      IF (A1.EQ.0.D0) GO TO 55
      DADR1=A1*(1./(R(1)-R(2))-1./(R(3)-R(1)))
      DADR2=A1*(1./(R(2)-R(3))-1./(R(1)-R(2)))
      DADR3=A1*(1./(R(3)-R(1))-1./(R(2)-R(3)))
      FA=-3.D0*ALF5*RS*RS*VA
      EX5=DEXP(-ALF5*RS**3)
      PAR=EX5*A1*(2.D0*A2+3.D0*A3*A1+4.D0*A4*A1*A1+5.D0*A5*A1*A1*A1) 
   55 DVA(1)=FA+PAR*DADR1  
      DVA(2)=FA+PAR*DADR2  
      DVA(3)=FA+PAR*DADR3  
      DB1DR1=(R(1)-R(2)-R(3))/A23-0.5D0*(S13/R(3)+S12/R(2))/R(1)**2  
      DB1DR2=(R(2)-R(1)-R(3))/A13-0.5D0*(S23/R(3)+S12/R(1))/R(2)**2  
      DB1DR3=(R(3)-R(1)-R(2))/A12-0.5D0*(S13/R(1)+S23/R(2))/R(3)**2  
      DB2DR1=-1./R(1)**2
      DB2DR2=-1./R(2)**2
      DB2DR3=-1./R(3)**2
      B3=(R(2)-R(1))**2+(R(3)-R(1))**2+(R(3)-R(2))**2
      DB3DR1=-2.*(R(2)+R(3)-R(1)-R(1)) 
      DB3DR2=-2.*(R(3)+R(1)-R(2)-R(2)) 
      DB3DR3=-2.*(R(1)+R(2)-R(3)-R(3)) 
      PARB=-BETA1*VB1+B1*B12*EXS1
      PARC=(B11+B12*RS)*EXS1
      DVA(1)=DVA(1)+PARB+PARC*DB1DR1
      DVA(2)=DVA(2)+PARB+PARC*DB1DR2
      DVA(3)=DVA(3)+PARB+PARC*DB1DR3
      PARB2=-2.*BETA2*RS*VB2
      PARB2K=B1*EXS2*(2.D0*B22+3.D0*B23*B1+4.D0*B24*B1*B1)
      DVA(1)=DVA(1)+PARB2+PARB2K*DB1DR1
      DVA(2)=DVA(2)+PARB2+PARB2K*DB1DR2
      DVA(3)=DVA(3)+PARB2+PARB2K*DB1DR3
      PAR32K=B1*(B31*EXS1+B32*EXS2*B1) 
      PAR31K=B2*(EXS1*B31+2.D0*B32*B1*EXS2)  
      PARB3=B2*(BETA1*B31*EXS1+2.D0*B32*BETA2*EXS2*RS*B1)*B1
      DVA(1)=DVA(1)+PAR32K*DB2DR1+PAR31K*DB1DR1-PARB3
      DVA(2)=DVA(2)+PAR32K*DB2DR2+PAR31K*DB1DR2-PARB3
      DVA(3)=DVA(3)+PAR32K*DB2DR3+PAR31K*DB1DR3-PARB3
      S14=B1*(B41*EXS1+B42*EXS2) 
      S34=B3*(B41*EXS1+B42*EXS2) 
      S134=-B1*B3*(BETA1*B41*EXS1+BETA2*B42*2.D0*EXS2*RS)
      DVA(1)=DVA(1)+DB3DR1*S14+DB1DR1*S34+S134
      DVA(2)=DVA(2)+DB3DR2*S14+DB1DR2*S34+S134
      DVA(3)=DVA(3)+DB3DR3*S14+DB1DR3*S34+S134
      S15=(B51+B52*RS*RS)*EXS3
      SSS5=2.*B1*B52*RS*EXS3-BETA3*VB5 
      DVA(1)=DVA(1)+S15*DB1DR1+SSS5
      DVA(2)=DVA(2)+S15*DB1DR2+SSS5
      DVA(3)=DVA(3)+S15*DB1DR3+SSS5
      DV(1) = (DQ(1)-((A(1)-A(2))*DA(1)-(A(3)-A(1))*DA(1))/
     +        (2.0*DSQRT(D)) + DVA(1)) 
      DV(2) = (DQ(2)-((A(2)-A(3))*DA(2)-(A(1)-A(2))*DA(2))/
     +        (2.0*DSQRT(D)) + DVA(2)) 
      DV(3) = (DQ(3)-((A(3)-A(1))*DA(3)-(A(2)-A(3))*DA(3))/
     +        (2.0*DSQRT(D)) + DVA(3)) 
      DV(1)=DV(1)*27.2107D0
      DV(2)=DV(2)*27.2107D0
      DV(3)=DV(3)*27.2107D0
      RETURN
      END

