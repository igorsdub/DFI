C****************************************************************
C	FORCED RESPONSE NETWORK PROGRAM
C****************************************************************
C	VERSION XXX 
C	WRITTEN/ARRANGED BY Z. Nevin Gerek, C. Atilgan
C
C	PLEASE REFER TO THE FOLLOWING PAPERS FOR MORE INFO:
C	
C	ANM:
C     Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin O, 
C     Bahar I, Biophys. J., 80:505-15, 2001. 
C
C	GNM: 
C	Bahar I, Atilgan AR, Erman B, Fold. & Des., 2:173-81, 1997
C
C****************************************************************
C	VARIABLES
C****************************************************************

C	NR: NUMBER OF RESIDUES 
C	CUTOFF: CUT-OFF RADIUS
C	EIGENCUT: CUT-OFF TO DECIDE ZERO EIGENVALUE(S)

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NR=171,GAMMA=100.0)
	REAL*8 X(NR),Y(NR),Z(NR),XF(NR),YF(NR),ZF(NR)
c    INTEGER*4 timeArray(3)    ! Holds the hour, minute, and second
	CHARACTER CNAM(NR)*3,CNAMF(NR)*3
	REAL*8 BETA(NR),HBETA(NR),conNum(NR),BETAF(NR)	
	REAL*8 BKBT(NR*3,NR*3)
	REAL*8 HESS(NR*3,NR*3),INVBKBT(NR*3,NR*3)
	REAL*8 W(NR*3),V(NR*3,NR*3)
	REAL*8 Wr(NR),Vr(NR,NR)
	REAL*8 Wrv(NR*3),Vrv(NR*3,NR*3)
	DIMENSION INDX(NR*3)
	REAL*8 INVCONT(NR,NR), sum2(NR), sum2r(NR*3)
	INTEGER RESNUM,RES3,ANUM(NR),RNUM(NR),ANUMF(NR),RNUMF(NR)
	INTEGER NUMZERO,STARTMOD,ENDMOD,totConMax,tIndex
	REAL*8 CUTOFFSQ,CUTOFF
	REAL*8 BX,BY,BZ,DIS2,EIGENCUT
	REAL*8 FLUCX2(NR),FLUCY2(NR),FLUCZ2(NR)
	REAL*8 FLUCX(NR),FLUCY(NR),FLUCZ(NR)
	REAL*8 FX2N(NR),FY2N(NR),FZ2N(NR)
	CHARACTER*4 ATNAME,ANAME(NR),ANAMEF(NR)
	CHARACTER*1 A1,CHA(NR),CHAF(NR)
	REAL*8 MEANB, MEANF,CORR,CORRF2,CORRB2
	REAL*8 delForce(NR*3),delRperb(NR*3),delR(NR*3,NR*3),delRf(NR*3),delrs(NR*3),delfs(NR*3)
	INTEGER nperturbRes,iperturbRes,fperturbRes
	REAL*8 XT(NR),YT(NR),ZT(NR),DIFF(NR),DTARGET(NR)
	REAL*8 XI(NR),YI(NR),ZI(NR)
	REAL*8 XTEMP(NR),YTEMP(NR),ZTEMP(NR),DIFBST(NR),rmsvalue,PR                
	REAL*8 delRperb_X(NR*3),delRperb_Y(NR*3),delRperb_Z(NR*3)
	REAL*8 delRperb_XY(NR*3),delRperb_YZ(NR*3),delRperb_XZ(NR*3)
	REAL*8 delRperb_XYZ(NR*3)
	REAL*8 PRCORR(NR)


	REAL*8 Xp(NR),Yp(NR),Zp(NR),Xn(NR),Yn(NR),Zn(NR)
	REAL*8 SQDIS2ratio_X(NR,NR), SQDIS2ratio_Y(NR,NR), SQDIS2ratio_Z(NR,NR)
	REAL*8 SQDIS2ratio_XY(NR,NR), SQDIS2ratio_XZ(NR,NR), SQDIS2ratio_YZ(NR,NR)
	
	REAL*8 SQDIS2ratio_XYZ(NR,NR), AVG(NR,NR), sumDIS2RATIO(NR)


C	DUMMIES

	INTEGER DINT,DINT2,ICA,IDUM
	REAL*8 DIFXX,DIFYY,DIFFZZ,DIST,DSUM,DSUM1,DSUM2
	CHARACTER DUMMY6*6,DUMMY3*3

	REAL*8 delRperbDist(NR),delRperbDistMat(NR,NR)


c *******************************************************************************
c        Declaration of variables for different directions respectively.
c *******************************************************************************

	REAL*8 XX(NR),YX(NR),ZX(NR)
	REAL*8 XT_X(NR),YT_X(NR),ZT_X(NR),DIFF_X(NR)
	REAL*8 XI_X(NR),YI_X(NR),ZI_X(NR)
	REAL*8 XTEMP_X(NR),YTEMP_X(NR),ZTEMP_X(NR),DIFBST_X(NR)

	REAL*8 XY(NR),YY(NR),ZY(NR)
	REAL*8 XT_Y(NR),YT_Y(NR),ZT_Y(NR),DIFF_Y(NR)
	REAL*8 XI_Y(NR),YI_Y(NR),ZI_Y(NR)
	REAL*8 XTEMP_Y(NR),YTEMP_Y(NR),ZTEMP_Y(NR),DIFBST_Y(NR)

	REAL*8 XZ(NR),YZ(NR),ZZ(NR)
	REAL*8 XT_Z(NR),YT_Z(NR),ZT_Z(NR),DIFF_Z(NR)
	REAL*8 XI_Z(NR),YI_Z(NR),ZI_Z(NR)
	REAL*8 XTEMP_Z(NR),YTEMP_Z(NR),ZTEMP_Z(NR),DIFBST_Z(NR)

	REAL*8 XXY(NR),YXY(NR),ZXY(NR)
	REAL*8 XT_XY(NR),YT_XY(NR),ZT_XY(NR),DIFF_XY(NR)
	REAL*8 XI_XY(NR),YI_XY(NR),ZI_XY(NR)
	REAL*8 XTEMP_XY(NR),YTEMP_XY(NR),ZTEMP_XY(NR),DIFBST_XY(NR)

	REAL*8 XXZ(NR),YXZ(NR),ZXZ(NR)
	REAL*8 XT_XZ(NR),YT_XZ(NR),ZT_XZ(NR),DIFF_XZ(NR)
	REAL*8 XI_XZ(NR),YI_XZ(NR),ZI_XZ(NR)
	REAL*8 XTEMP_XZ(NR),YTEMP_XZ(NR),ZTEMP_XZ(NR),DIFBST_XZ(NR)

	REAL*8 XYZ(NR),YYZ(NR),ZYZ(NR)
	REAL*8 XT_YZ(NR),YT_YZ(NR),ZT_YZ(NR),DIFF_YZ(NR)
	REAL*8 XI_YZ(NR),YI_YZ(NR),ZI_YZ(NR)
	REAL*8 XTEMP_YZ(NR),YTEMP_YZ(NR),ZTEMP_YZ(NR),DIFBST_YZ(NR)

	REAL*8 XXYZ(NR),YXYZ(NR),ZXYZ(NR)
	REAL*8 XT_XYZ(NR),YT_XYZ(NR),ZT_XYZ(NR),DIFF_XYZ(NR)
	REAL*8 XI_XYZ(NR),YI_XYZ(NR),ZI_XYZ(NR)
	REAL*8 XTEMP_XYZ(NR),YTEMP_XYZ(NR),ZTEMP_XYZ(NR),DIFBST_XYZ(NR)


	REAL*8 delRperbDist_X(NR),delRperbDistMat_X(NR,NR),delRperbDistMatlog_X(NR,NR)
	REAL*8 delRperbDist_Y(NR),delRperbDistMat_Y(NR,NR),delRperbDistMatlog_Y(NR,NR)
	REAL*8 delRperbDist_Z(NR),delRperbDistMat_Z(NR,NR),delRperbDistMatlog_Z(NR,NR)
	REAL*8 delRperbDist_XY(NR),delRperbDistMat_XY(NR,NR),delRperbDistMatlog_XY(NR,NR)
	REAL*8 delRperbDist_YZ(NR),delRperbDistMat_YZ(NR,NR),delRperbDistMatlog_YZ(NR,NR)
	REAL*8 delRperbDist_XZ(NR),delRperbDistMat_XZ(NR,NR),delRperbDistMatlog_XZ(NR,NR)
	REAL*8 delRperbDist_XYZ(NR),delRperbDistMat_XYZ(NR,NR),delRperbDistMatlog_XYZ(NR,NR)
	REAL*8 delRperbDistMat_Avg(NR,NR)
	REAL*8 Dmag(NR,NR),DmagRow(NR),ds(NR)
	
	REAL*8 sum1ixlog(NR),sum1iylog(NR),sum1izlog(NR),sum1ixylog(NR),sum1iyzlog(NR),sum1ixzlog(NR),sum1ixyzlog(NR)
	REAL*8 sum1jxlog(NR),sum1jylog(NR),sum1jzlog(NR),sum1jxylog(NR),sum1jyzlog(NR),sum1jxzlog(NR),sum1jxyzlog(NR)
	REAL*8 sum1ix(NR),sum1iy(NR),sum1iz(NR),sum1ixy(NR),sum1iyz(NR),sum1ixz(NR),sum1ixyz(NR)
	REAL*8 sum1jx(NR),sum1jy(NR),sum1jz(NR),sum1jxy(NR),sum1jyz(NR),sum1jxz(NR),sum1jxyz(NR)
	REAL*8 AvglogRow(NR),AvglogCol(NR),logAvgnormCol(NR),AvgCol(NR),AvgnormCol(NR),AvgColA
	REAL*8 AvgRow(NR),AvgnormRow(NR)
	REAL*8 delRperbDistMat_XNorm(NR,NR),sum2ix(NR),sum2ixlog(NR),plog5X(NR)
	REAL*8 delRperbDistMat_YNorm(NR,NR),sum2iy(NR),sum2iylog(NR),plog5Y(NR)
	REAL*8 delRperbDistMat_ZNorm(NR,NR),sum2iz(NR),sum2izlog(NR),plog5Z(NR)
	REAL*8 delRperbDistMat_XYNorm(NR,NR),sum2ixy(NR),sum2ixylog(NR),plog5XY(NR)
	REAL*8 delRperbDistMat_XZNorm(NR,NR),sum2ixz(NR),sum2ixzlog(NR),plog5XZ(NR)
	REAL*8 delRperbDistMat_YZNorm(NR,NR),sum2iyz(NR),sum2iyzlog(NR),plog5YZ(NR)
	REAL*8 delRperbDistMat_XYZNorm(NR,NR),sum2ixyz(NR),sum2ixyzlog(NR),plog5XYZ(NR)
	REAL*8 corrcoeff(NR),test1(NR),u12(NR),u22(NR)
	REAL*8 delRperbMat_XYZNorm(NR*3,NR*3),sum5(NR*3),delRperbMat_XYZ(NR*3,NR*3)
	REAL*8 tSUM1(NR),tSUM3(NR)
	REAL*8 INVBKBTRes(3,3),vRes(3,3),wRes(3),traceInvBKBT(NR),delRperbMat_XYZRes(3,3)
	REAL*8 delRperb_X_X(NR,NR), delRperb_Y_X(NR,NR), delRperb_Z_X(NR,NR)
	REAL*8 delRperb_X_Y(NR,NR), delRperb_Y_Y(NR,NR), delRperb_Z_Y(NR,NR)
	REAL*8 delRperb_X_Z(NR,NR), delRperb_Y_Z(NR,NR), delRperb_Z_Z(NR,NR)
	REAL*8 delRperb_X_XY(NR,NR),delRperb_Y_XY(NR,NR),delRperb_Z_XY(NR,NR)
	REAL*8 delRperb_X_XZ(NR,NR),delRperb_Y_XZ(NR,NR),delRperb_Z_XZ(NR,NR)
	REAL*8 delRperb_X_YZ(NR,NR),delRperb_Y_YZ(NR,NR),delRperb_Z_YZ(NR,NR)
	REAL*8 delRperb_X_XYZ(NR,NR),delRperb_Y_XYZ(NR,NR),delRperb_Z_XYZ(NR,NR)
	
	REAL*8 sum_X_X(NR),sum_Y_X(NR),sum_Z_X(NR)
	REAL*8 sum_X_Y(NR),sum_Y_Y(NR),sum_Z_Y(NR)
	REAL*8 sum_X_Z(NR),sum_Y_Z(NR),sum_Z_Z(NR)
	REAL*8 sum_X_XY(NR),sum_Y_XY(NR),sum_Z_XY(NR)
	REAL*8 sum_X_XZ(NR),sum_Y_XZ(NR),sum_Z_XZ(NR)
	REAL*8 sum_X_YZ(NR),sum_Y_YZ(NR),sum_Z_YZ(NR)
	REAL*8 sum_X_XYZ(NR),sum_Y_XYZ(NR),sum_Z_XYZ(NR)
	REAL*8 distS2_X(NR),distS2_Y(NR),distS2_Z(NR),distS2_XY(NR),distS2_XZ(NR),distS2_YZ(NR),distS2_XYZ(NR)
	REAL*8 MaxCol(NR),MaxnormCol(NR)
	REAL*8 MaxRow(NR),MaxnormRow(NR),MaxdistS2(NR),AvgdistS2(NR)


C******************************************************************
C	PARAMETERS
C******************************************************************

	RESNUM=NR
	RES3=NR*3

C	radius threshold to decide if two residues/nucleotides are 
C	connected
	OPEN(99,FILE='input.parms')
	READ(99,*) iperturbRes,fperturbRes,perbegin,perend,alpha
	EIGENCUT=1E-6

C******************************************************************
C	FILES
C******************************************************************

C	THESE ARE THE INPUT FILES
C	files can be obtained from Brookhaven Protein Data Bank 
C	http://www.rcsb.org/pdb
C	the 1st is the structure to be perturbed, the 2nd is the target 
	OPEN(50,FILE='uf.pdb')
	OPEN(51,FILE='uf.pdb')

c *******output files

c	open(20,file='uf-hessianPF.txt')

	OPEN(60,FILE='init_centers.pdb')
	OPEN(59,FILE='target_centers.pdb')
	open(61,file='uf.delRinitial')
	OPEN(66,FILE='eigenvalues.txt')
	open(69,file='uf_RESULTS.txt')
	OPEN(71,FILE='uf_diffbest.txt')
	OPEN(72,FILE='uf_dev_corr.txt')


	OPEN (95,file='20sloweigenx.txt')
	OPEN (96,file='20sloweigeny.txt')
	OPEN (97,file='20sloweigenz.txt')
	OPEN (98,file='20sloweigenr.txt')

	
	open(88,file='uf-S2-Avg.dat')
	open(89,file='uf-S2-Max.dat')
	open(90,file='uf-S1-Avg.dat')
	open(91,file='uf-S1-Max.dat')

	open(15,file='uf-ResponseVector-Magnitude-Avg.dat')
	open(16,file='uf-ResponseVector-Magnitude-Max.dat')
	

C******************************************************************
C	READ COORDINATES, AND B-FACTORS - INITIAL FILE
C******************************************************************
	
310	READ(50,'(A6)') DUMMY6
	IF(DUMMY6.NE.'ATOM  ') GOTO 310
	BACKSPACE(50)

	ICA=1
320	READ(50,'(A6)') DUMMY6
	IF(DUMMY6.NE.'ATOM  ') THEN
	 IF(DUMMY6.EQ.'END   ') THEN
	  GOTO 330
	 ELSE
	  GOTO 320
	 ENDIF
	ENDIF
	BACKSPACE(50)
      READ(50,55) DUMMY6,D1INT,ATNAME,DUMMY3,A1,D2INT,XXX,YYY,ZZZ,R1,BBB
	
	

	IF(ATNAME.EQ.' CA ') THEN
	 ANUM(ICA)=D1INT
	 ANAME(ICA)=ATNAME
	 CNAM(ICA)=DUMMY3
	 CHA(ICA)=A1
	 RNUM(ICA)=D2INT
	 X(ICA)=XXX
	 Y(ICA)=YYY
	 Z(ICA)=ZZZ
	 BETA(ICA)=BBB
	 ICA=ICA+1
	ENDIF
	GOTO 320
	
55	FORMAT(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,2F6.2)

	WRITE(*,*) RESNUM

330	IF(RESNUM.NE.(ICA-1)) THEN 
	 WRITE(*,*) 'THERE IS A PROBLEM WITH THE NUMBER OF RESIDUES!'
	 WRITE(*,*) 'GIVEN RESNUM=',RESNUM,'CALCULATED RESNUM=',ICA-1
	 GOTO 666
	ENDIF

	WRITE(60,*)"INITIAL COORDINATES"
	DO 10 I=1,RESNUM
10    write(60,55)"ATOM  ",I,ANAME(I),CNAM(I),"A",I,X(I),Y(I),Z(I),
     +RNUM(I),BETA(I)
	write(60,56)'END'

56	FORMAT(A3)


C******************************************************************
C	READ COORDINATES, AND B-FACTORS - TARGET FILE
C******************************************************************

311	READ(51,'(A6)') DUMMY6
	IF(DUMMY6.NE.'ATOM  ') GOTO 311
	BACKSPACE(51)

	ICA=1
321	READ(51,'(A6)') DUMMY6
	IF(DUMMY6.NE.'ATOM  ') THEN
	 IF(DUMMY6.EQ.'END   ') THEN
	  GOTO 331
	 ELSE
	  GOTO 321
	 ENDIF
	ENDIF
	BACKSPACE(51)
      READ(51,55) DUMMY6,D1INT,ATNAME,DUMMY3,A1,D2INT,XXX,YYY,ZZZ,R1,BBB

C     INTERCHANGE THE FOLLOWING COMMENT LINES IF NETWORK IS TO BE CENTERED ON CB
	IF(ATNAME.EQ.' CA ') THEN
c	IF(ATNAME.EQ.' CB '.OR.(ATNAME.EQ.' CA '.AND.DUMMY3.EQ.'GLY'))THEN
	 ANUMF(ICA)=D1INT
	 ANAMEF(ICA)=ATNAME
	 CNAMF(ICA)=DUMMY3
	 CHAF(ICA)=A1
	 RNUMF(ICA)=D2INT
	 XF(ICA)=XXX
	 YF(ICA)=YYY
	 ZF(ICA)=ZZZ
	 BETAF(ICA)=BBB
	 ICA=ICA+1
	ENDIF
	GOTO 321
	
331	IF(RESNUM.NE.(ICA-1)) THEN 
	 WRITE(*,*) 'THERE IS A PROBLEM WITH THE NUMBER OF RESIDUES!'
	 WRITE(*,*) 'GIVEN RESNUM=',RESNUM,'CALCULATED RESNUM=',ICA-1
	 GOTO 666
	ENDIF

	WRITE(59,*)"TARGET COORDINATES"
	DO 11 I=1,RESNUM
11    write(59,55)"ATOM  ",I,ANAMEF(I),CNAMF(I),"A",I,XF(I),YF(I),ZF(I),
     +RNUMF(I),BETAF(I)
	write(59,56)'END'


C	SUPERIMPOSE THE TWO STRUCTURES; DTARGET HOLDS THE TARGET DIFFERENCES
C	NOTE THAT IN subroutine impose, IF SPECIFIC IMMOBILE DOMAINS ARE TO 
C	BE SUPERPOSED THAN THOSE RESIDUES MUST BE ENTERED IN THE SUBROUTINE:
c	superposition will best-fit on atoms for wfit=1.

	DO I=1,RESNUM
	XI(I) = X(I)
	YI(I) = Y(I)
	ZI(I) = Z(I)
	ENDDO


	call impose(NR,XI,YI,ZI,NR,XF,YF,ZF,rmsvalue)
	

	ISAY=0
	sumdtarget =0.0
	DO I=1,RESNUM
		DX=XF(I)-XI(I)
		DY=YF(I)-YI(I)
		DZ=ZF(I)-ZI(I)
		DTARGET(I)=SQRT(DX*DX+DY*DY+DZ*DZ)	
		sumdtarget = sumdtarget + DTARGET(I)
	ENDDO

	do I=1,RESNUM	
		write(61,*)I, (XF(I)-XI(I)),(YF(I)-YI(I)),(ZF(I)-ZI(I))
	enddo



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C	THIS IS THE OUTER MAIN LOOP, WHERE THE STRUCTURE CHOSEN FROM
C	THE PREVIOUS STEP IS PERTURBED RESIDUE-BY-RESIDUE
C	CONTINUE PERTURBING UNTIL THE TARGET DIFFERENCES ARE REACHED
C	OR AN UPPER LIMIT FOR THE NUMBER OF ITERATIONS (ITERF)
C
C	FOR A SINGLE STEP, LET ITERF=0
C
C	DELF IS THE MAXIMUM SIZE OF THE FORCE APPLIED 
C	one has the option of perturbing in a non-random manner; see DELF loop below
C	PRTARGET IS THE TARGETED PEARSON CORRELATION
C	RMSTARGET IS THE TARGETED RMSD BETWEEN THE PERTURBED STRUCTURE AND THE EXPERIMENTAL.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	ITERF=0
	DELF=0.1
	PRTARGET=0.999
	RMSTARGET=0.1
c	subroutine ran1 uses needs a negative integer seed for initialization:
	idum=-67126379

500	CONTINUE

	ISAY=ISAY+1

	write(72,*)"STEP ",ISAY
	write(66,*)"STEP ",ISAY

C	PR0 will store the highest Pearson correlation in each iteration
C	RMS0 will store the lowest RMSD in each iteration

	PR0=0.
	RMS0=1000000.


C******************************************************************
C	determine the number of connected interactions
C******************************************************************
      DO J=1,RESNUM-1
	  conNum(J) = 0
	  DO K=J+1,RESNUM
		BX=X(J)-X(K)
		BY=Y(J)-Y(K)
		BZ=Z(J)-Z(K)
		DIS2=BX*BX+BY*BY+BZ*BZ
        IF(J.NE.K) THEN
c        IF(J.NE.K.AND.DIS2.LE.CUTOFFSQ) THEN
			 conNum(J) =conNum(J) +1
		ENDIF			
	ENDDO
	ENDDO
	
	sum=0
	DO J=1, RESNUM
		sum =sum+ conNum(J)
	ENDDO
	
	totConMax = sum
	
c	 WRITE(*,*) 'NO. OF MEMBERS:',totConMax
	 

C******************************************************************
C	INITIALIZATION OF HESSIAN MATRIX
C******************************************************************

	DO 20 I=1,RES3
	DO 20 J=1,RES3
20	 HESS(I,J)=0.

C******************************************************************
C	CREATION OF HESSIAN MATRIX
C******************************************************************

	DO J=1,RESNUM
	DO K=1,RESNUM
	 BX=X(J)-X(K)
	 BY=Y(J)-Y(K)
	 BZ=Z(J)-Z(K)
         DIS2=BX*BX+BY*BY+BZ*BZ

       IF(J.NE.K) THEN
C	  FIRST: CREATION OF Hii
	  HESS(3*J-2,3*J-2)=HESS(3*J-2,3*J-2)+(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BX/DIS2
	  HESS(3*J-1,3*J-1)=HESS(3*J-1,3*J-1)+(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BY/DIS2
	  HESS(3*J,3*J)=HESS(3*J,3*J)        +(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BZ*BZ/DIS2

	  HESS(3*J-2,3*J-1)=HESS(3*J-2,3*J-1)+(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BY/DIS2
	  HESS(3*J-2,3*J)=HESS(3*J-2,3*J)    +(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BZ/DIS2
	  HESS(3*J-1,3*J-2)=HESS(3*J-1,3*J-2)+(GAMMA /DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BX/DIS2
	  HESS(3*J-1,3*J)=HESS(3*J-1,3*J)    +(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BZ/DIS2
	  HESS(3*J,3*J-2)=HESS(3*J,3*J-2)    +(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BZ/DIS2
	  HESS(3*J,3*J-1)=HESS(3*J,3*J-1)    +(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BZ/DIS2
       
C	  SECOND: CREATION OF  Hij
	  HESS(3*J-2,3*K-2)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BX/DIS2
	  HESS(3*J-1,3*K-1)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BY/DIS2
	  HESS(3*J,3*K)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BZ*BZ/DIS2

	  HESS(3*J-2,3*K-1)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BY/DIS2
	  HESS(3*J-2,3*K)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BZ/DIS2
	  HESS(3*J-1,3*K-2)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BX/DIS2
	  HESS(3*J-1,3*K)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BZ/DIS2
	  HESS(3*J,3*K-2)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BX*BZ/DIS2
	  HESS(3*J,3*K-1)=	-(GAMMA/DIS2)*(GAMMA/DIS2)*(GAMMA/DIS2)*BY*BZ/DIS2 
       ENDIF			
	ENDDO
	ENDDO	 


C******************************************************************
C	SINGULAR VALUE DECOMPOSITION TO GET RID OF ZERO EIGENVALUES
C******************************************************************

	CALL SVDCMP(HESS,RES3,RES3,RES3,RES3,W,V)

C******************************************************************
C	PUTTING THE EIGENVALUES IN ASCENDING ORDER
C******************************************************************

	CALL INDEXX(RES3,W,INDX)

	DO 700 I=1,RES3
	 WRITE(66,*) I,W(INDX(I))
700	CONTINUE

C******************************************************************
C	SAVING THE EIGENVECTORS
C******************************************************************
	NEIG=RES3
	LN=RES3

C******************************************************************
C	CALCULATING ZERO EIGENVALUES
C******************************************************************
	NUMZERO=0
	DO  K=1,RES3
	 IF(W(K).LE.EIGENCUT) THEN
		NUMZERO=NUMZERO+1
	 ENDIF
	ENDDO
	 WRITE(*,*) 'NUMBER OF ZERO EIGENVALUES',NUMZERO

C	NUMZERO should be equal to 6 due the our connnectivity definition
	IF(NUMZERO.NE.6)WRITE(*,*) 'IRREGULAR ZERO EIGENVALUES!!!',NUMZERO

C******************************************************************
C	Calculating inverse of Hessian Matrix
C******************************************************************

	
	DO I=1,RES3
	 DO J=1,RES3
	  INVBKBT(I,J)=0.
	  DO  K=1,RES3
	   IF(W(K).GT.EIGENCUT) THEN
	    INVBKBT(I,J)=INVBKBT(I,J)+V(I,K)*HESS(J,K)/W(K)
         ENDIF
	  ENDDO
       ENDDO
	ENDDO


	JJJ=0
	DO I=1,RES3
	 IF(INVBKBT(I,I).LT.0) JJJ=JJJ+1
      ENDDO
	IF(JJJ.GT.0)WRITE(*,*) 'NEGATIVE ELEMENTS IN THE DIAGONAL:',JJJ



C******************************************************************
C	MODESHAPES - Twenty Slowest 
C******************************************************************
 
		k1 = NUMZERO + 1
		k2 =  NUMZERO + 20
		do 170 i=1,RES3, 3
			ii = (i-1)/3 + 1
			write(95,21273)ii,(v(i,INDX(k)),k=k1,k2,1)
			write(96,21273)ii,(v(i+1,INDX(k)),k=k1,k2,1)
			write(97,21273)ii,(v(i+2,INDX(k)),k=k1,k2,1)
			write(98,21273)ii,(v(i,INDX(k))*v(i,INDX(k))+v(i+1,INDX(k))*v(i+1,INDX(k))+v(i+2,INDX(k))*v(i+2,INDX(k)),k=k1,k2,1)
170		continue
21273	format(i7,22f11.6)		


C******************************************************************
C	Define force vector - sequential perturbation over all residues
C	THIS IS THE MAIN INNER LOOP
C	EACH RESIDUE IS PERTURBED AND THE STRUCTURE WITH THE BEST 
C	(i) DIFFERENCE FROM THE TARGET STRUCTURE (as measured by
C	Pearson CORRELATION, (ii) rmsd FROM THE TARGET STRUCTURE, or
C	(iii) A COMBINATION OF THE TWO
C	IS SELECTED AS THE NEXT CANDIDATE TO BE PERTUBED
C******************************************************************


c********************************************************************
c********************************************************************

c    Perturbation through X-direction
	
c********************************************************************
c********************************************************************

	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO I=1,RESNUM
		Xp(I) = 0.
 		Yp(I) = 0.
		Zp(I) = 0.
		Xn(I) = 0.
		Yn(I) = 0.
		Zn(I) = 0.
	ENDDO
	
	DO 1000 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
		  delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
		DO i=3*nperturbRes-2,3*nperturbRes-2
c		   delForce(i) = DELF*(ran1(idum)*2-1)
		   delForce(i) = 1.0
		ENDDO
		
C	Multiply force vector with inverse of Hessian = DeltaR
	 
	DO I=1,RES3
	 delRperb_X(I)= 0.0
	 DO J=1,RES3
	    delRperb_X(I)= delRperb_X(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO

     
	DO I=1,RESNUM
	  XI_X(I) = XI(I)
      YI_X(I) = YI(I)
      ZI_X(I) = ZI(I)	
	ENDDO

c calculate new coordinates based on any ratio (add and subtract)
	DO I=perbegin,perend
	   Xp(I) = XI_X(I) + alpha * delRperb_X(3*I-2)
	   Yp(I) = YI_X(I) + alpha * delRperb_X(3*I-1)
	   Zp(I) = ZI_X(I) + alpha * delRperb_X(3*I)
	ENDDO

	DO I=perbegin,perend
	   Xn(I) = XI_X(I) - alpha * delRperb_X(3*I-2)
	   Yn(I) = YI_X(I) - alpha * delRperb_X(3*I-1)
	   Zn(I) = ZI_X(I) - alpha * delRperb_X(3*I)
	ENDDO

c write into output files
c	if(nperturbRes.lt.10) then
c		write(25,411) nperturbRes
c		write(26,411) nperturbRes
c		write(75,411) nperturbRes
c411		FORMAT('PerturbRes',I1)
c	elseif(nperturbRes.lt.100) then
c		write(25,412) nperturbRes
c		write(26,412) nperturbRes
c		write(75,412) nperturbRes
c412		FORMAT('PerturbRes',I2)
c	else
c		write(25,413) nperturbRes
c		write(26,413) nperturbRes
c		write(75,413) nperturbRes
c413		FORMAT('PerturbRes',I3)
c	endif

c	do I=1,RESNUM
c		write(25,*)I, Xp(I),Yp(I),Zp(I)
c		write(26,*)I, Xn(I),Yn(I),Zn(I)
c	enddo

	do I=1,RESNUM	
	  XT_X(I) = 0.
	  YT_X(I) = 0.
	  ZT_X(I) = 0.
      enddo

	do I=1,RESNUM	
		XT_X(I) = XI_X(I) +  delRperb_X(3*I-2) 
		YT_X(I) = YI_X(I) +  delRperb_X(3*I-1) 
		ZT_X(I) = ZI_X(I) +  delRperb_X(3*I) 
	enddo

C	Superimpose the new structure to the initial,

	call impose(NR,XI_X,YI_X,ZI_X,NR,XT_X,YT_X,ZT_X,rmsvalue)


C	COMPUTE THE DIFFERENCES

	DO I=1,RESNUM
		DX_X=XT_X(I)-XI_X(I)
		DY_X=YT_X(I)-YI_X(I)
		DZ_X=ZT_X(I)-ZI_X(I)
		DIFF_X(I)=SQRT(DX_X*DX_X+DY_X*DY_X+DZ_X*DZ_X)
	ENDDO

	call impose(NR,XF,YF,ZF,NR,XT_X,YT_X,ZT_X,rmsvalue)
	CALL pearsn (DIFF_X,DTARGET,NR,PR)

	WRITE(72,*)'X direction'
c	WRITE(*,*)perturbRes,rmsvalue,PR
	WRITE(72,*)nperturbRes,rmsvalue,PR

	CALL INDEXX(RESNUM,PRCORR,INDX)


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
c	IF(PR.GT.PR0)THEN
c	IF(PR.GT.PR0.AND.rmsvalue.LT.RMS0)THEN
C	IF(rmsvalue.LT.RMS0)THEN
		PR0=PR
		RMS0=rmsvalue
		ICHOSEN=nperturbRes
		DO I=1,RESNUM
			XTEMP_X(I) = XT_X(I) 
			YTEMP_X(I) = YT_X(I) 
			ZTEMP_X(I) = ZT_X(I) 
			DIFBST_X(I) = DIFF_X(I)
		ENDDO
c	ENDIF
	
	do I=perbegin,perend
	   ratioXDiff_X = alpha * (XTEMP_X(I)-XI_X(I))
	   ratioYDiff_X = alpha * (YTEMP_X(I)-YI_X(I))
	   ratioZDiff_X = alpha * (ZTEMP_X(I)-ZI_X(I))
	   SQDIS2ratio_X(nperturbRes,I) = sqrt ( ratioXDiff_X*ratioXDiff_X + ratioYDiff_X * ratioYDiff_X + ratioZDiff_X * ratioZDiff_X )
c		write(75,*)I, ratioXDiff_X,ratioYDiff_X,ratioZDiff_X,SQDIS2ratio_X(nperturbRes,I)
	enddo

c ------------------

	   DO i=1,RESNUM
		   delRperbDist_X(i) = sqrt( delRperb_X(3*I-2) *  delRperb_X(3*I-2) +  delRperb_X(3*I-1) * delRperb_X(3*I-1) +  delRperb_X(3*I) *  delRperb_X(3*I) )
	   ENDDO
 
	   DO i=nperturbRes,nperturbRes
	    DO j=1,RESNUM
			delRperbDistMat_X(i,j) = delRperbDist_X(j)
			delRperb_X_X(i,j) = delRperb_X(3*J-2) 
			delRperb_Y_X(i,j) = delRperb_X(3*J-1) 
			delRperb_Z_X(i,j) = delRperb_X(3*J) 
		ENDDO
	   ENDDO


c---------------------


1000	enddo


c---  delRperb_X_X,delRperb_Y_X,delRperb_Z_X wrt row

	   do i=1,resnum
	   sum_X_X(i) = 0.0
	   sum_Y_X(i) = 0.0
	   sum_Z_X(i) = 0.0
c	   write(18,*) 	'Perturbed Residue', i
	   do j=1,resnum	  
c			write(18,*) i,j,delRperb_X_XYZ(i,j)
	   		sum_X_X(i) = sum_X_X(i) + delRperb_X_X(i,j)
	   		sum_Y_X(i) = sum_Y_X(i) + delRperb_Y_X(i,j)
	   		sum_Z_X(i) = sum_Z_X(i) + delRperb_Z_X(i,j)
	   enddo
			sum_X_X(i) = sum_X_X(i)/resnum
			sum_Y_X(i) = sum_Y_X(i)/resnum
			sum_Z_X(i) = sum_Z_X(i)/resnum
	   enddo
	   

c	   do i=1,resnum
c			write(14,*) sum_X(i),sum_Y(i),sum_Z(i)
c	   enddo


	   sumNormX=0.0
	   do i=1,resnum
		distS2_X(i) = sqrt(sum_X_X(i)*sum_X_X(i)+sum_Y_X(i)*sum_Y_X(i)+sum_Z_X(i)*sum_Z_X(i))
		sumNormX = sumNormX + distS2_X(i)
	   enddo
	   do i=1,resnum
		 distS2_X(i)= distS2_X(i)/sumNormX
	   enddo




c--- normalize matrix wrt row sum


	   sum3=0.0
	   do i=1,resnum
	   sum2(i) = 0.0
	   do j=1,resnum
	   		sum2(i) = sum2(i) + delRperbDistMat_X(i,j)
	   enddo
	   sum3 = sum3 + sum2(i)
	   enddo
	   
	   do i=1,resnum
c	   		write(*,*) i, sum2(i)
	   		do j=1,resnum
	   			delRperbDistMat_XNorm(i,j) = delRperbDistMat_X(i,j) /sum2(i)
	   		enddo
	   enddo


c column sum
	   DO i=1,RESNUM
	   sum1jx(i) = 0.0
		DO j=1,RESNUM
			sum1jx(i) = sum1jx(i) + delRperbDistMat_X(j,i)
		ENDDO
c		write(*,*) i,sum1jx(i)
	   ENDDO


c row sum
	   DO i=1,RESNUM
	   sum1ix(i) = 0.0
		DO j=1,RESNUM
			sum1ix(i) = sum1ix(i) + delRperbDistMat_X(i,j)
		ENDDO
c		write(*,*) i,sum1jx(i)
	   ENDDO


	
	WRITE(69,*)'X direction'
	WRITE(*,511)ISAY,RMS0,PR0,ICHOSEN
	WRITE(69,511)ISAY,RMS0,PR0,ICHOSEN
511	FORMAT('ITER:',I4,3X,'RMS:',F7.3,3X,'CORR:',F6.3,3X,'RES:',I4)

	RMS0X = RMS0
	DO I=1,RESNUM
		XX(I) = XTEMP_X(I) 
		YX(I) = YTEMP_X(I) 
		ZX(I) = ZTEMP_X(I) 
	ENDDO

	WRITE(71,*)'X direction'
	WRITE(71,*)'STEP',ISAY
	do I=1,RESNUM	
      write(71,55)"ATOM  ",I," CA ",CNAM(I),"A",I,XX(I),YX(I),ZX(I)
     +	  ,RNUM(I),BETA(I)
	WRITE(71,*)I,DIFBST_X(I)
	enddo
	write(71,56)'END'


c********************************************************************
c********************************************************************

c                            Y-direction
	
c********************************************************************
c********************************************************************

	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO I=1,RESNUM
		Xp(I) = 0.
 		Yp(I) = 0.
		Zp(I) = 0.
		Xn(I) = 0.
		Yn(I) = 0.
		Zn(I) = 0.
	ENDDO
    	
	DO 1001 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO

C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
		DO i=3*nperturbRes-1,3*nperturbRes-1
c			delForce(i) = DELF*(ran1(idum)*2-1)
			delForce(i) = 1.0
		ENDDO

C	Multiply force vector with inverse of Hessian = DeltaR
	 
	DO I=1,RES3
	 delRperb_Y(I)= 0.0
	 DO J=1,RES3
	    delRperb_Y(I)= delRperb_Y(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO

	DO I=1,RESNUM
	    XI_Y(I) = XI(I)
        YI_Y(I) = YI(I)
        ZI_Y(I) = ZI(I)	
	ENDDO
   
c calculate new coordinates based on any ratio
	DO I=perbegin,perend
	   Xp(I) = XI_Y(I) + alpha * delRperb_Y(3*I-2)
	   Yp(I) = YI_Y(I) + alpha * delRperb_Y(3*I-1)
	   Zp(I) = ZI_Y(I) + alpha * delRperb_Y(3*I)
	ENDDO

	DO I=perbegin,perend
	   Xn(I) = XI_Y(I) - alpha * delRperb_Y(3*I-2)
	   Yn(I) = YI_Y(I) - alpha * delRperb_Y(3*I-1)
	   Zn(I) = ZI_Y(I) - alpha * delRperb_Y(3*I)
	ENDDO

c write into output files
c	if(nperturbRes.lt.10) then
c		write(29,411) nperturbRes
c		write(30,411) nperturbRes
c		write(76,411) nperturbRes
c	elseif(nperturbRes.lt.100) then
c		write(29,412) nperturbRes
c		write(30,412) nperturbRes
c		write(76,412) nperturbRes
c	else
c		write(29,413) nperturbRes
c		write(30,413) nperturbRes
c		write(76,413) nperturbRes
c	endif

c	do I=1,RESNUM
c		write(29,*)I, Xp(I),Yp(I),Zp(I)
c		write(30,*)I, Xn(I),Yn(I),Zn(I)
c	enddo

	DO I=1,RESNUM
	   XT_Y(I) = 0.
       YT_Y(I) = 0.
       ZT_Y(I) = 0.	
	ENDDO

	do I=1,RESNUM	
		XT_Y(I) = XI_Y(I) +  delRperb_Y(3*I-2) 
		YT_Y(I) = YI_Y(I) +  delRperb_Y(3*I-1) 
		ZT_Y(I) = ZI_Y(I) +  delRperb_Y(3*I) 
	enddo

C	Superimpose the new structure to the initial,

	call impose(NR,XI_Y,YI_Y,ZI_Y,NR,XT_Y,YT_Y,ZT_Y,rmsvalue)


C	COMPUTE THE DIFFERENCES

	DO I=1,RESNUM
		DX_Y=XT_Y(I)-XI_Y(I)
		DY_Y=YT_Y(I)-YI_Y(I)
		DZ_Y=ZT_Y(I)-ZI_Y(I)
		DIFF_Y(I)=SQRT(DX_Y*DX_Y+DY_Y*DY_Y+DZ_Y*DZ_Y)
	ENDDO

	call impose(NR,XF,YF,ZF,NR,XT_Y,YT_Y,ZT_Y,rmsvalue)
	CALL pearsn (DIFF_Y,DTARGET,NR,PR)

	WRITE(72,*)'Y diRection'
c	WRITE(*,*)perturbRes,rmsvalue,PR
	WRITE(72,*)nperturbRes,rmsvalue,PR

	CALL INDEXX(RESNUM,PRCORR,INDX)


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
c	IF(PR.GT.PR0)THEN
c	IF(PR.GT.PR0.AND.rmsvalue.LT.RMS0)THEN
C	IF(rmsvalue.LT.RMS0)THEN
		PR0=PR
		RMS0=rmsvalue
		ICHOSEN=nperturbRes
		DO I=1,RESNUM
			XTEMP_Y(I) = XT_Y(I) 
			YTEMP_Y(I) = YT_Y(I) 
			ZTEMP_Y(I) = ZT_Y(I) 
			DIFBST_Y(I) = DIFF_Y(I)
		ENDDO
c	ENDIF

	do I=perbegin,perend
	   ratioXDiff_Y = alpha * (XTEMP_Y(I)-XI_Y(I))
	   ratioYDiff_Y = alpha * (YTEMP_Y(I)-YI_Y(I))
	   ratioZDiff_Y = alpha * (ZTEMP_Y(I)-ZI_Y(I))
	   SQDIS2ratio_Y(nperturbRes,I) = sqrt ( ratioXDiff_Y*ratioXDiff_Y + ratioYDiff_Y * ratioYDiff_Y + ratioZDiff_Y * ratioZDiff_Y )
c		write(76,*)I, ratioXDiff_Y,ratioYDiff_Y,ratioZDiff_Y,SQDIS2ratio_Y(nperturbRes,I)
	enddo
c ------------------

	   DO i=1,RESNUM
		   delRperbDist_Y(i) = sqrt( delRperb_Y(3*I-2) *  delRperb_Y(3*I-2) +  delRperb_Y(3*I-1) * delRperb_Y(3*I-1) +  delRperb_Y(3*I) *  delRperb_Y(3*I) )
	   ENDDO
	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_Y(i,j) = delRperbDist_Y(j) 
			delRperb_X_Y(i,j) = delRperb_Y(3*J-2) 
			delRperb_Y_Y(i,j) = delRperb_Y(3*J-1) 
			delRperb_Z_Y(i,j) = delRperb_Y(3*J) 
		ENDDO
	   ENDDO


c---------------------


1001	enddo


c---  delRperb_X_X,delRperb_Y_X,delRperb_Z_X wrt row

	   do i=1,resnum
	   sum_X_Y(i) = 0.0
	   sum_Y_Y(i) = 0.0
	   sum_Z_Y(i) = 0.0
c	   write(18,*) 	'Perturbed Residue', i
	   do j=1,resnum	  
c			write(18,*) i,j,delRperb_X_XYZ(i,j)
	   		sum_X_Y(i) = sum_X_Y(i) + delRperb_X_Y(i,j)
	   		sum_Y_Y(i) = sum_Y_Y(i) + delRperb_Y_Y(i,j)
	   		sum_Z_Y(i) = sum_Z_Y(i) + delRperb_Z_Y(i,j)
	   enddo
			sum_X_Y(i) = sum_X_Y(i)/resnum
			sum_Y_Y(i) = sum_Y_Y(i)/resnum
			sum_Z_Y(i) = sum_Z_Y(i)/resnum
	   enddo
	   

c	   do i=1,resnum
c			write(14,*) sum_X(i),sum_Y(i),sum_Z(i)
c	   enddo


	   sumNormY=0.0
	   do i=1,resnum
		distS2_Y(i) = sqrt(sum_X_Y(i)*sum_X_Y(i)+sum_Y_Y(i)*sum_Y_Y(i)+sum_Z_Y(i)*sum_Z_Y(i))
		sumNormY = sumNormY + distS2_Y(i)
	   enddo
	   do i=1,resnum
		 distS2_Y(i)= distS2_Y(i)/sumNormY
	   enddo



c--- normalize matrix based on row

	   sum3=0.0
	   do i=1,resnum
	   sum2(i) = 0.0
	   do j=1,resnum
	   		sum2(i) = sum2(i) + delRperbDistMat_Y(i,j)
	   enddo
	   sum3 = sum3 + sum2(i)
	   enddo
	   
	   do i=1,resnum
	   	do j=1,resnum
	   		delRperbDistMat_YNorm(i,j) = delRperbDistMat_Y(i,j) /sum2(i)
	   	enddo
	   enddo

c column sum
	   DO i=1,RESNUM
	   sum1jy(i) = 0.0
		DO j=1,RESNUM
			sum1jy(i) = sum1jy(i) + delRperbDistMat_Y(j,i)
		ENDDO
	   ENDDO
c row sum
	   DO i=1,RESNUM
	   sum1iy(i) = 0.0
		DO j=1,RESNUM
			sum1iy(i) = sum1iy(i) + delRperbDistMat_Y(i,j)
		ENDDO
c		write(*,*) i,sum1iy(i)
	   ENDDO

	
	WRITE(69,*)'Y direction'
	WRITE(*,511)ISAY,RMS0,PR0,ICHOSEN
	WRITE(69,511)ISAY,RMS0,PR0,ICHOSEN

	RMS0X = RMS0
	DO I=1,RESNUM
		XY(I) = XTEMP_Y(I) 
		YY(I) = YTEMP_Y(I) 
		ZY(I) = ZTEMP_Y(I) 
	ENDDO

	WRITE(71,*)'Y direction'
	WRITE(71,*)'STEP',ISAY
	do I=1,RESNUM	
      write(71,55)"ATOM  ",I," CA ",CNAM(I),"A",I,XY(I),YY(I),ZY(I)
     +	  ,RNUM(I),BETA(I)
	WRITE(71,*)I,DIFBST_Y(I)
	enddo
	write(71,56)'END'


c********************************************************************
c********************************************************************

c  Perturbation through Z-direction
	
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO I=1,RESNUM
		Xp(I) = 0.
 		Yp(I) = 0.
		Zp(I) = 0.
		Xn(I) = 0.
		Yn(I) = 0.
		Zn(I) = 0.
	ENDDO
	DO 1002 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes,3*nperturbRes
c			delForce(i) = DELF*(ran1(idum)*2-1)
			delForce(i) = 1.0
		ENDDO
		
C	Multiply force vector with inverse of Hessian = DeltaR
	 
	DO I=1,RES3
	 delRperb_Z(I)= 0.0
	 DO J=1,RES3
	    delRperb_Z(I)= delRperb_Z(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO

	DO I=1,RESNUM
	    XI_Z(I) = XI(I)
        YI_Z(I) = YI(I)
        ZI_Z(I) = ZI(I)	
	ENDDO
 
c calculate new coordinates based on any ratio
	DO I=perbegin,perend
	   Xp(I) = XI_Z(I) + alpha * delRperb_Z(3*I-2)
	   Yp(I) = YI_Z(I) + alpha * delRperb_Z(3*I-1)
	   Zp(I) = ZI_Z(I) + alpha * delRperb_Z(3*I)
	ENDDO

	DO I=perbegin,perend
	   Xn(I) = XI_Z(I) - alpha * delRperb_Z(3*I-2)
	   Yn(I) = YI_Z(I) - alpha * delRperb_Z(3*I-1)
	   Zn(I) = ZI_Z(I) - alpha * delRperb_Z(3*I)
	ENDDO

c	if(nperturbRes.lt.10) then
c		write(33,411) nperturbRes
c		write(34,411) nperturbRes
c		write(77,411) nperturbRes
c	elseif(nperturbRes.lt.100) then
c		write(33,412) nperturbRes
c		write(34,412) nperturbRes
c		write(77,412) nperturbRes
c	else
c		write(33,413) nperturbRes
c		write(34,413) nperturbRes
c		write(77,413) nperturbRes
c	endif

c	do I=1,RESNUM
c		write(33,*)I, Xp(I),Yp(I),Zp(I)
c		write(34,*)I, Xn(I),Yn(I),Zn(I)
c	enddo
    

	DO I=1,RESNUM
 		  XT_Z(I) = 0.
          YT_Z(I) = 0.
          ZT_Z(I) = 0.	
	ENDDO

	do I=1,RESNUM	
	   XT_Z(I) = XI_Z(I) +  delRperb_Z(3*I-2) 
	   YT_Z(I) = YI_Z(I) +  delRperb_Z(3*I-1) 
	   ZT_Z(I) = ZI_Z(I) +  delRperb_Z(3*I) 
	enddo

C	Superimpose the new structure to the initial,

	call impose(NR,XI_Z,YI_Z,ZI_Z,NR,XT_Z,YT_Z,ZT_Z,rmsvalue)


C	COMPUTE THE DIFFERENCES

	DO I=1,RESNUM
	   DX_Z=XT_Z(I)-XI_Z(I)
	   DY_Z=YT_Z(I)-YI_Z(I)
	   DZ_Z=ZT_Z(I)-ZI_Z(I)
	   DIFF_Z(I)=SQRT(DX_Z*DX_Z+DY_Z*DY_Z+DZ_Z*DZ_Z)
	ENDDO

	call impose(NR,XF,YF,ZF,NR,XT_Z,YT_Z,ZT_Z,rmsvalue)
	CALL pearsn (DIFF_Z,DTARGET,NR,PR)

	WRITE(72,*)'Z diRection'
c	WRITE(*,*)perturbRes,rmsvalue,PR
	WRITE(72,*)nperturbRes,rmsvalue,PR

	CALL INDEXX(RESNUM,PRCORR,INDX)


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
c	IF(PR.GT.PR0)THEN
c	IF(PR.GT.PR0.AND.rmsvalue.LT.RMS0)THEN
C	IF(rmsvalue.LT.RMS0)THEN
		PR0=PR
		RMS0=rmsvalue
		ICHOSEN=nperturbRes
		DO I=1,RESNUM
			XTEMP_Z(I) = XT_Z(I) 
			YTEMP_Z(I) = YT_Z(I) 
			ZTEMP_Z(I) = ZT_Z(I) 
			DIFBST_Z(I) = DIFF_Z(I)
		ENDDO
c	ENDIF

	do I=perbegin,perend
	   ratioXDiff_Z = alpha * (XTEMP_Z(I)-XI_Z(I))
	   ratioYDiff_Z = alpha * (YTEMP_Z(I)-YI_Z(I))
	   ratioZDiff_Z = alpha * (ZTEMP_Z(I)-ZI_Z(I))
	   SQDIS2ratio_Z(nperturbRes,I) = sqrt ( ratioXDiff_Z*ratioXDiff_Z + ratioYDiff_Z * ratioYDiff_Z + ratioZDiff_Z * ratioZDiff_Z )
c		write(77,*)I, ratioXDiff_Z,ratioYDiff_Z,ratioZDiff_Z,SQDIS2ratio_Z(nperturbRes,I)
	enddo

c ------------------


	   DO i=1,RESNUM
		   delRperbDist_Z(i) = sqrt( delRperb_Z(3*I-2) *  delRperb_Z(3*I-2) +  delRperb_Z(3*I-1) * delRperb_Z(3*I-1) +  delRperb_Z(3*I) *  delRperb_Z(3*I) )
	   ENDDO
	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_Z(i,j) = delRperbDist_Z(j) 
			delRperb_X_Z(i,j) = delRperb_Z(3*J-2) 
			delRperb_Y_Z(i,j) = delRperb_Z(3*J-1) 
			delRperb_Z_Z(i,j) = delRperb_Z(3*J) 
		ENDDO
	   ENDDO


c---------------------


1002	enddo


c---  delRperb_X_X,delRperb_Y_X,delRperb_Z_X wrt row

	   do i=1,resnum
	   sum_X_Z(i) = 0.0
	   sum_Y_Z(i) = 0.0
	   sum_Z_Z(i) = 0.0
c	   write(18,*) 	'Perturbed Residue', i
	   do j=1,resnum	  
c			write(18,*) i,j,delRperb_X_XYZ(i,j)
	   		sum_X_Z(i) = sum_X_Z(i) + delRperb_X_Z(i,j)
	   		sum_Y_Z(i) = sum_Y_Z(i) + delRperb_Y_Z(i,j)
	   		sum_Z_Z(i) = sum_Z_Z(i) + delRperb_Z_Z(i,j)
	   enddo
			sum_X_Z(i) = sum_X_Z(i)/resnum
			sum_Y_Z(i) = sum_Y_Z(i)/resnum
			sum_Z_Z(i) = sum_Z_Z(i)/resnum
	   enddo
	   

c	   do i=1,resnum
c			write(14,*) sum_X(i),sum_Y(i),sum_Z(i)
c	   enddo


	   sumNormZ=0.0
	   do i=1,resnum
		distS2_Z(i) = sqrt(sum_X_Z(i)*sum_X_Z(i)+sum_Y_Z(i)*sum_Y_Z(i)+sum_Z_Z(i)*sum_Z_Z(i))
		sumNormZ = sumNormZ + distS2_Z(i)
	   enddo
	   do i=1,resnum
		 distS2_Z(i)= distS2_Z(i)/sumNormZ
	   enddo



c--- normalize matrix based on row

	   sum3=0.0
	   do i=1,resnum
	   sum2(i) = 0.0
	   do j=1,resnum
	   		sum2(i) = sum2(i) + delRperbDistMat_Z(i,j)
	   enddo
	   sum3 = sum3 + sum2(i)
	   enddo
	   
	   do i=1,resnum
	   	do j=1,resnum
	   		delRperbDistMat_ZNorm(i,j) = delRperbDistMat_Z(i,j) /sum2(i)
	   	enddo
	   enddo

c column sum
	   DO i=1,RESNUM
	   sum1jz(i) = 0.0
		DO j=1,RESNUM
			sum1jz(i) = sum1jz(i) + delRperbDistMat_Z(j,i)
		ENDDO
	   ENDDO
c row sum
	   DO i=1,RESNUM
	   sum1iz(i) = 0.0
		DO j=1,RESNUM
			sum1iz(i) = sum1iz(i) + delRperbDistMat_Z(i,j)
		ENDDO
	   ENDDO

	   
	WRITE(69,*)'Z direction'
	WRITE(*,511)ISAY,RMS0,PR0,ICHOSEN
	WRITE(69,511)ISAY,RMS0,PR0,ICHOSEN

	RMS0X = RMS0
	DO I=1,RESNUM
		XZ(I) = XTEMP_Z(I) 
		YZ(I) = YTEMP_Z(I) 
		ZZ(I) = ZTEMP_Z(I) 
	ENDDO

	WRITE(71,*)'Z direction'
	WRITE(71,*)'STEP',ISAY
	do I=1,RESNUM	
      write(71,55)"ATOM  ",I," CA ",CNAM(I),"A",I,XZ(I),YZ(I),ZZ(I)
     +	  ,RNUM(I),BETA(I)
	WRITE(71,*)I,DIFBST_Z(I)
	enddo
	write(71,56)'END'


c********************************************************************
c********************************************************************

c   Perturbation through XY-direction
	
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO I=1,RESNUM
		Xp(I) = 0.
 		Yp(I) = 0.
		Zp(I) = 0.
		Xn(I) = 0.
		Yn(I) = 0.
		Zn(I) = 0.
	ENDDO
    
	DO 1003 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-2,3*nperturbRes-1
c			delForce(i) = DELF*(ran1(idum)*2-1)
			delForce(i) = 1.0/sqrt(2.0)
		ENDDO
		
C	Multiply force vector with inverse of Hessian = DeltaR
	 
	DO I=1,RES3
	 delRperb_XY(I)= 0.0
	 DO J=1,RES3
	    delRperb_XY(I)= delRperb_XY(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
	    XI_XY(I) = XI(I)
        YI_XY(I) = YI(I)
        ZI_XY(I) = ZI(I)	
	ENDDO

c calculate new coordinates based on any ratio
	DO I=perbegin,perend
	   Xp(I) = XI_XY(I) + alpha * delRperb_XY(3*I-2)
	   Yp(I) = YI_XY(I) + alpha * delRperb_XY(3*I-1)
	   Zp(I) = ZI_XY(I) + alpha * delRperb_XY(3*I)
	ENDDO

	DO I=perbegin,perend
	   Xn(I) = XI_XY(I) - alpha * delRperb_XY(3*I-2)
	   Yn(I) = YI_XY(I) - alpha * delRperb_XY(3*I-1)
	   Zn(I) = ZI_XY(I) - alpha * delRperb_XY(3*I)
	ENDDO

c write into output files
c	if(nperturbRes.lt.10) then
c		write(37,411) nperturbRes
c		write(38,411) nperturbRes
c		write(78,411) nperturbRes
c	elseif(nperturbRes.lt.100) then
c		write(37,412) nperturbRes
c		write(38,412) nperturbRes
c		write(78,412) nperturbRes
c	else
c		write(37,413) nperturbRes
c		write(38,413) nperturbRes
c		write(78,413) nperturbRes
c
c	endif

c	do I=1,RESNUM
c		write(37,*)I, Xp(I),Yp(I),Zp(I)
c		write(38,*)I, Xn(I),Yn(I),Zn(I)
c	enddo

	DO I=1,RESNUM
	      XT_XY(I) = 0.
          YT_XY(I) = 0.
          ZT_XY(I) = 0.	
	ENDDO

	do I=1,RESNUM	
	  XT_XY(I) = XI_XY(I) +  delRperb_XY(3*I-2) 
	  YT_XY(I) = YI_XY(I) +  delRperb_XY(3*I-1) 
	  ZT_XY(I) = ZI_XY(I) +  delRperb_XY(3*I) 
	enddo

C	Superimpose the new structure to the initial,

	call impose(NR,XI_XY,YI_XY,ZI_XY,NR,XT_XY,YT_XY,ZT_XY,rmsvalue)


C	COMPUTE THE DIFFERENCES

	DO I=1,RESNUM
		DX_XY=XT_XY(I)-XI_XY(I)
		DY_XY=YT_XY(I)-YI_XY(I)
		DZ_XY=ZT_XY(I)-ZI_XY(I)
		DIFF_XY(I)=SQRT(DX_XY*DX_XY+DY_XY*DY_XY+DZ_XY*DZ_XY)
	ENDDO

	call impose(NR,XF,YF,ZF,NR,XT_XY,YT_XY,ZT_XY,rmsvalue)
	CALL pearsn (DIFF_XY,DTARGET,NR,PR)

	WRITE(72,*)'XY diRection'
c	WRITE(*,*)perturbRes,rmsvalue,PR
	WRITE(72,*)nperturbRes,rmsvalue,PR

	CALL INDEXX(RESNUM,PRCORR,INDX)


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
c	IF(PR.GT.PR0)THEN
c	IF(PR.GT.PR0.AND.rmsvalue.LT.RMS0)THEN
C	IF(rmsvalue.LT.RMS0)THEN
		PR0=PR
		RMS0=rmsvalue
		ICHOSEN=nperturbRes
		DO I=1,RESNUM
			XTEMP_XY(I) = XT_XY(I) 
			YTEMP_XY(I) = YT_XY(I) 
			ZTEMP_XY(I) = ZT_XY(I) 
			DIFBST_XY(I) = DIFF_XY(I)
		ENDDO
c	ENDIF

	do I=perbegin,perend
	   ratioXDiff_XY = alpha * (XTEMP_XY(I)-XI_XY(I))
	   ratioYDiff_XY = alpha * (XTEMP_XY(I)-YI_XY(I))
	   ratioZDiff_XY = alpha * (XTEMP_XY(I)-ZI_XY(I))
	   SQDIS2ratio_XY(nperturbRes,I) = sqrt ( ratioXDiff_XY*ratioXDiff_XY + ratioYDiff_XY * ratioYDiff_XY + ratioZDiff_XY * ratioZDiff_XY )
c		write(78,*)I, ratioXDiff_XY,ratioYDiff_XY,ratioZDiff_XY,SQDIS2ratio_XY(nperturbRes,I)
	enddo

c ------------------

	   DO i=1,RESNUM
		   delRperbDist_XY(i) = sqrt( delRperb_XY(3*I-2) *  delRperb_XY(3*I-2) +  delRperb_XY(3*I-1) * delRperb_XY(3*I-1) +  delRperb_XY(3*I) *  delRperb_XY(3*I) )
	   ENDDO
	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_XY(i,j) = delRperbDist_XY(j) 
			delRperb_X_XY(i,j) = delRperb_XY(3*J-2) 
			delRperb_Y_XY(i,j) = delRperb_XY(3*J-1) 
			delRperb_Z_XY(i,j) = delRperb_XY(3*J) 
		ENDDO
	   ENDDO


c---------------------


1003	enddo

c---  delRperb_X_X,delRperb_Y_X,delRperb_Z_X wrt row

	   do i=1,resnum
	   sum_X_XY(i) = 0.0
	   sum_Y_XY(i) = 0.0
	   sum_Z_XY(i) = 0.0
c	   write(18,*) 	'Perturbed Residue', i
	   do j=1,resnum	  
c			write(18,*) i,j,delRperb_X_XYZ(i,j)
	   		sum_X_XY(i) = sum_X_XY(i) + delRperb_X_XY(i,j)
	   		sum_Y_XY(i) = sum_Y_XY(i) + delRperb_Y_XY(i,j)
	   		sum_Z_XY(i) = sum_Z_XY(i) + delRperb_Z_XY(i,j)
	   enddo
			sum_X_XY(i) = sum_X_XY(i)/resnum
			sum_Y_XY(i) = sum_Y_XY(i)/resnum
			sum_Z_XY(i) = sum_Z_XY(i)/resnum
	   enddo
	   

c	   do i=1,resnum
c			write(14,*) sum_X(i),sum_Y(i),sum_Z(i)
c	   enddo


	   sumNormXY=0.0
	   do i=1,resnum
		distS2_XY(i) = sqrt(sum_X_XY(i)*sum_X_XY(i)+sum_Y_XY(i)*sum_Y_XY(i)+sum_Z_XY(i)*sum_Z_XY(i))
		sumNormXY = sumNormXY + distS2_XY(i)
	   enddo
	   do i=1,resnum
		 distS2_XY(i)= distS2_XY(i)/sumNormXY
	   enddo



c--- normalize matrix based on row

	   sum3=0.0
	   do i=1,resnum
	   sum2(i) = 0.0
	   do j=1,resnum
	   		sum2(i) = sum2(i) + delRperbDistMat_XY(i,j)
	   enddo
	   sum3 = sum3 + sum2(i)
	   enddo
	   
	   do i=1,resnum
	   	do j=1,resnum
	   		delRperbDistMat_XYNorm(i,j) = delRperbDistMat_XY(i,j) /sum2(i)
	   	enddo
	   enddo

c column sum
	   DO i=1,RESNUM
	   sum1jxy(i) = 0.0
		DO j=1,RESNUM
			sum1jxy(i) = sum1jxy(i) + delRperbDistMat_XY(j,i)
		ENDDO
	   ENDDO
	   
c row sum
	   DO i=1,RESNUM
	   sum1ixy(i) = 0.0
		DO j=1,RESNUM
			sum1ixy(i) = sum1ixy(i) + delRperbDistMat_XY(i,j)
		ENDDO
	   ENDDO

	   	
	WRITE(69,*)'XY direction'
	WRITE(*,511)ISAY,RMS0,PR0,ICHOSEN
	WRITE(69,511)ISAY,RMS0,PR0,ICHOSEN

	RMS0X = RMS0
	DO I=1,RESNUM
		XXY(I) = XTEMP_XY(I) 
		YXY(I) = YTEMP_XY(I) 
		ZXY(I) = ZTEMP_XY(I) 
	ENDDO

	WRITE(71,*)'XY direction'
	WRITE(71,*)'STEP',ISAY
	do I=1,RESNUM	
      write(71,55)"ATOM  ",I," CA ",CNAM(I),"A",I,XXY(I),YXY(I),ZXY(I)
     +	  ,RNUM(I),BETA(I)
	WRITE(71,*)I,DIFBST_XY(I)
	enddo
	write(71,56)'END'


c********************************************************************
c********************************************************************

c Perturbation through XZ-direction
	
c********************************************************************
c********************************************************************

	DO i=1,NR
	  PRCORR(i) = 0.
	ENDDO

	DO I=1,RESNUM
		Xp(I) = 0.
 		Yp(I) = 0.
		Zp(I) = 0.
		Xn(I) = 0.
		Yn(I) = 0.
		Zn(I) = 0.
	ENDDO
    
	DO 1004 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-2,3*nperturbRes,2
c			delForce(i) = DELF*(ran1(idum)*2-1)
			delForce(i) = 1.0/sqrt(2.0)
		ENDDO
		

C	Multiply force vector with inverse of Hessian = DeltaR


	DO I=1,RES3
	 delRperb_XZ(I)= 0.0
	 DO J=1,RES3
	    delRperb_XZ(I)= delRperb_XZ(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
	    XI_XZ(I) = XI(I)
        YI_XZ(I) = YI(I)
        ZI_XZ(I) = ZI(I)	
	ENDDO

c calculate new coordinates based on any ratio
	DO I=perbegin,perend
	   Xp(I) = XI_XZ(I) + alpha * delRperb_XZ(3*I-2)
	   Yp(I) = YI_XZ(I) + alpha * delRperb_XZ(3*I-1)
	   Zp(I) = ZI_XZ(I) + alpha * delRperb_XZ(3*I)
	ENDDO

	DO I=perbegin,perend
	   Xn(I) = XI_XZ(I) - alpha * delRperb_XZ(3*I-2)
	   Yn(I) = YI_XZ(I) - alpha * delRperb_XZ(3*I-1)
	   Zn(I) = ZI_XZ(I) - alpha * delRperb_XZ(3*I)
	ENDDO
c write into output files
c	if(nperturbRes.lt.10) then
c		write(41,411) nperturbRes
c		write(42,411) nperturbRes
c		write(79,411) nperturbRes
c	elseif(nperturbRes.lt.100) then
c		write(41,412) nperturbRes
c		write(42,412) nperturbRes
c		write(79,412) nperturbRes
c	else
c		write(41,413) nperturbRes
c		write(42,413) nperturbRes
c		write(79,413) nperturbRes
c	endif

c	do I=1,RESNUM
c		write(41,*)I, Xp(I),Yp(I),Zp(I)
c		write(42,*)I, Xn(I),Yn(I),Zn(I)
c	enddo


	DO I=1,RESNUM
		 XT_XZ(I) = 0.
         YT_XZ(I) = 0.
         ZT_XZ(I) = 0.	
	ENDDO

	do I=1,RESNUM	
		XT_XZ(I) = XI_XZ(I) +  delRperb_XZ(3*I-2) 
		YT_XZ(I) = YI_XZ(I) +  delRperb_XZ(3*I-1) 
		ZT_XZ(I) = ZI_XZ(I) +  delRperb_XZ(3*I) 
	enddo

C	Superimpose the new structure to the initial,

	call impose(NR,XI_XZ,YI_XZ,ZI_XZ,NR,XT_XZ,YT_XZ,ZT_XZ,rmsvalue)


C	COMPUTE THE DIFFERENCES

	DO I=1,RESNUM
		DX_XZ=XT_XZ(I)-XI_XZ(I)
		DY_XZ=YT_XZ(I)-YI_XZ(I)
		DZ_XZ=ZT_XZ(I)-ZI_XZ(I)
		DIFF_XZ(I)=SQRT(DX_XZ*DX_XZ+DY_XZ*DY_XZ+DZ_XZ*DZ_XZ)
	ENDDO

	call impose(NR,XF,YF,ZF,NR,XT_XZ,YT_XZ,ZT_XZ,rmsvalue)
	CALL pearsn (DIFF_XZ,DTARGET,NR,PR)

	WRITE(72,*)'XZ diRection'
c	WRITE(*,*)perturbRes,rmsvalue,PR
	WRITE(72,*)nperturbRes,rmsvalue,PR

	CALL INDEXX(RESNUM,PRCORR,INDX)


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
c	IF(PR.GT.PR0)THEN
c	IF(PR.GT.PR0.AND.rmsvalue.LT.RMS0)THEN
C	IF(rmsvalue.LT.RMS0)THEN
		PR0=PR
		RMS0=rmsvalue
		ICHOSEN=nperturbRes
		DO I=1,RESNUM
			XTEMP_XZ(I) = XT_XZ(I) 
			YTEMP_XZ(I) = YT_XZ(I) 
			ZTEMP_XZ(I) = ZT_XZ(I) 
			DIFBST_XZ(I) = DIFF_XZ(I)
		ENDDO
c	ENDIF

	do I=perbegin,perend
	   ratioXDiff_XZ = alpha * (XTEMP_XZ(I)-XI_XZ(I))
	   ratioYDiff_XZ = alpha * (XTEMP_XZ(I)-YI_XZ(I))
	   ratioZDiff_XZ = alpha * (XTEMP_XZ(I)-ZI_XZ(I))
	   SQDIS2ratio_XZ(nperturbRes,I) = sqrt ( ratioXDiff_XZ*ratioXDiff_XZ + ratioYDiff_XZ * ratioYDiff_XZ + ratioZDiff_XZ * ratioZDiff_XZ )
c		write(79,*)I, ratioXDiff_XZ,ratioYDiff_XZ,ratioZDiff_XZ,SQDIS2ratio_XZ(nperturbRes,I)
	enddo

c ------------------

	   DO i=1,RESNUM
		   delRperbDist_XZ(i) = sqrt( delRperb_XZ(3*I-2) *  delRperb_XZ(3*I-2) +  delRperb_XZ(3*I-1) * delRperb_XZ(3*I-1) +  delRperb_XZ(3*I) *  delRperb_XZ(3*I) )
	   ENDDO
	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_XZ(i,j) = delRperbDist_XZ(j) 
			delRperb_X_XZ(i,j) = delRperb_XZ(3*J-2) 
			delRperb_Y_XZ(i,j) = delRperb_XZ(3*J-1) 
			delRperb_Z_XZ(i,j) = delRperb_XZ(3*J) 
		ENDDO
	   ENDDO


c---------------------


1004	enddo
	
c---  delRperb_X_X,delRperb_Y_X,delRperb_Z_X wrt row

	   do i=1,resnum
	   sum_X_XZ(i) = 0.0
	   sum_Y_XZ(i) = 0.0
	   sum_Z_XZ(i) = 0.0
c	   write(18,*) 	'Perturbed Residue', i
	   do j=1,resnum	  
c			write(18,*) i,j,delRperb_X_XYZ(i,j)
	   		sum_X_XZ(i) = sum_X_XZ(i) + delRperb_X_XZ(i,j)
	   		sum_Y_XZ(i) = sum_Y_XZ(i) + delRperb_Y_XZ(i,j)
	   		sum_Z_XZ(i) = sum_Z_XZ(i) + delRperb_Z_XZ(i,j)
	   enddo
			sum_X_XZ(i) = sum_X_XZ(i)/resnum
			sum_Y_XZ(i) = sum_Y_XZ(i)/resnum
			sum_Z_XZ(i) = sum_Z_XZ(i)/resnum
	   enddo
	   

c	   do i=1,resnum
c			write(14,*) sum_X(i),sum_Y(i),sum_Z(i)
c	   enddo


	   sumNormXZ=0.0
	   do i=1,resnum
		distS2_XZ(i) = sqrt(sum_X_XZ(i)*sum_X_XZ(i)+sum_Y_XZ(i)*sum_Y_XZ(i)+sum_Z_XZ(i)*sum_Z_XZ(i))
		sumNormXZ = sumNormXZ + distS2_XZ(i)
	   enddo
	   do i=1,resnum
		 distS2_XZ(i)= distS2_XZ(i)/sumNormXZ
	   enddo



c--- normalize matrix based on row

	   sum3=0.0
	   do i=1,resnum
	   sum2(i) = 0.0
	   do j=1,resnum
	   		sum2(i) = sum2(i) + delRperbDistMat_XZ(i,j)
	   enddo
	   sum3 = sum3 + sum2(i)
	   enddo
	   
	   do i=1,resnum
	   	do j=1,resnum
	   		delRperbDistMat_XZNorm(i,j) = delRperbDistMat_XZ(i,j) /sum2(i)
	   	enddo
	   enddo

c  column sum
	   DO i=1,RESNUM
	   sum1jxz(i) = 0.0
		DO j=1,RESNUM
			sum1jxz(i) = sum1jxz(i) + delRperbDistMat_XZ(j,i)
		ENDDO
	   ENDDO

c  row sum
	   DO i=1,RESNUM
	   sum1ixz(i) = 0.0
		DO j=1,RESNUM
			sum1ixz(i) = sum1ixz(i) + delRperbDistMat_XZ(i,j)
		ENDDO
	   ENDDO

	   
	WRITE(69,*)'XZ direction'
	WRITE(*,511)ISAY,RMS0,PR0,ICHOSEN
	WRITE(69,511)ISAY,RMS0,PR0,ICHOSEN

	RMS0X = RMS0
	DO I=1,RESNUM
		XXZ(I) = XTEMP_XZ(I) 
		YXZ(I) = YTEMP_XZ(I) 
		ZXZ(I) = ZTEMP_XZ(I) 
	ENDDO

	WRITE(71,*)'XZ direction'
	WRITE(71,*)'STEP',ISAY
	do I=1,RESNUM	
      write(71,55)"ATOM  ",I," CA ",CNAM(I),"A",I,XXZ(I),YXZ(I),ZXZ(I)
     +	  ,RNUM(I),BETA(I)
	WRITE(71,*)I,DIFBST_XZ(I)
	enddo
	write(71,56)'END'


c********************************************************************
c********************************************************************

c Perturbation through YZ-direction
	
c********************************************************************
c********************************************************************

	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO I=1,RESNUM
		Xp(I) = 0.
 		Yp(I) = 0.
		Zp(I) = 0.
		Xn(I) = 0.
		Yn(I) = 0.
		Zn(I) = 0.
	ENDDO
    
	DO 1005 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-1,3*nperturbRes
c			delForce(i) = DELF*(ran1(idum)*2-1)
			delForce(i) = 1.0/sqrt(2.0)
		ENDDO
		
C	Multiply force vector with inverse of BBT = DeltaR
	 
	DO I=1,RES3
	 delRperb_YZ(I)= 0.0
	 DO J=1,RES3
	    delRperb_YZ(I)= delRperb_YZ(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
	    XI_YZ(I) = XI(I)
        YI_YZ(I) = YI(I)
        ZI_YZ(I) = ZI(I)	
	ENDDO

c calculate new coordinates based on any ratio
	DO I=perbegin,perend
	   Xp(I) = XI_YZ(I) + alpha * delRperb_YZ(3*I-2)
	   Yp(I) = YI_YZ(I) + alpha * delRperb_YZ(3*I-1)
	   Zp(I) = ZI_YZ(I) + alpha * delRperb_YZ(3*I)
	ENDDO

	DO I=perbegin,perend
	   Xn(I) = XI_YZ(I) - alpha * delRperb_YZ(3*I-2)
	   Yn(I) = YI_YZ(I) - alpha * delRperb_YZ(3*I-1)
	   Zn(I) = ZI_YZ(I) - alpha * delRperb_YZ(3*I)
	ENDDO

c	if(nperturbRes.lt.10) then
c		write(45,411) nperturbRes
c		write(46,411) nperturbRes
c		write(80,411) nperturbRes
c	elseif(nperturbRes.lt.100) then
c		write(45,412) nperturbRes
c		write(46,412) nperturbRes
c		write(80,412) nperturbRes
c	else
c		write(45,413) nperturbRes
c		write(46,413) nperturbRes
c		write(80,413) nperturbRes
c	endif
c
c	do I=1,RESNUM
c		write(45,*)I, Xp(I),Yp(I),Zp(I)
c		write(46,*)I, Xn(I),Yn(I),Zn(I)
c	enddo


	DO I=1,RESNUM
		 XT_YZ(I) = 0.
         YT_YZ(I) = 0.
         ZT_YZ(I) = 0.	
	ENDDO

	do I=1,RESNUM	
		XT_YZ(I) = XI_YZ(I) +  delRperb_YZ(3*I-2) 
		YT_YZ(I) = YI_YZ(I) +  delRperb_YZ(3*I-1) 
		ZT_YZ(I) = ZI_YZ(I) +  delRperb_YZ(3*I) 
	enddo

C	Superimpose the new structure to the initial,

	call impose(NR,XI_YZ,YI_YZ,ZI_YZ,NR,XT_YZ,YT_YZ,ZT_YZ,rmsvalue)


C	COMPUTE THE DIFFERENCES

	DO I=1,RESNUM
		DX_YZ=XT_YZ(I)-XI_YZ(I)
		DY_YZ=YT_YZ(I)-YI_YZ(I)
		DZ_YZ=ZT_YZ(I)-ZI_YZ(I)
		DIFF_YZ(I)=SQRT(DX_YZ*DX_YZ+DY_YZ*DY_YZ+DZ_YZ*DZ_YZ)
	ENDDO

	call impose(NR,XF,YF,ZF,NR,XT_YZ,YT_YZ,ZT_YZ,rmsvalue)
	CALL pearsn (DIFF_YZ,DTARGET,NR,PR)

	WRITE(72,*)'YZ diRection'
c	WRITE(*,*)perturbRes,rmsvalue,PR
	WRITE(72,*)nperturbRes,rmsvalue,PR

	CALL INDEXX(RESNUM,PRCORR,INDX)


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
c	IF(PR.GT.PR0)THEN
c	IF(PR.GT.PR0.AND.rmsvalue.LT.RMS0)THEN
C	IF(rmsvalue.LT.RMS0)THEN
		PR0=PR
		RMS0=rmsvalue
		ICHOSEN=nperturbRes
		DO I=1,RESNUM
			XTEMP_YZ(I) = XT_YZ(I) 
			YTEMP_YZ(I) = YT_YZ(I) 
			ZTEMP_YZ(I) = ZT_YZ(I) 
			DIFBST_YZ(I) = DIFF_YZ(I)
		ENDDO
c	ENDIF

	do I=perbegin,perend
	   ratioXDiff_YZ = alpha * (XTEMP_YZ(I)-XI_YZ(I))
	   ratioYDiff_YZ = alpha * (XTEMP_YZ(I)-YI_YZ(I))
	   ratioZDiff_YZ = alpha * (XTEMP_YZ(I)-ZI_YZ(I))
	   SQDIS2ratio_YZ(nperturbRes,I) = sqrt ( ratioXDiff_YZ*ratioXDiff_YZ + ratioYDiff_YZ * ratioYDiff_YZ + ratioZDiff_YZ * ratioZDiff_YZ )
c		write(80,*)I, ratioXDiff_YZ,ratioYDiff_YZ,ratioZDiff_YZ,SQDIS2ratio_YZ(nperturbRes,I)
	enddo

c ------------------

	   DO i=1,RESNUM
		   delRperbDist_YZ(i) = sqrt( delRperb_YZ(3*I-2) *  delRperb_YZ(3*I-2) +  delRperb_YZ(3*I-1) * delRperb_YZ(3*I-1) +  delRperb_YZ(3*I) *  delRperb_YZ(3*I) )
	   ENDDO
	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_YZ(i,j) = delRperbDist_YZ(j) 
			delRperb_X_YZ(i,j) = delRperb_YZ(3*J-2) 
			delRperb_Y_YZ(i,j) = delRperb_YZ(3*J-1) 
			delRperb_Z_YZ(i,j) = delRperb_YZ(3*J) 
		ENDDO
	   ENDDO


c---------------------


1005	enddo

c---  delRperb_X_X,delRperb_Y_X,delRperb_Z_X wrt row

	   do i=1,resnum
	   sum_X_YZ(i) = 0.0
	   sum_Y_YZ(i) = 0.0
	   sum_Z_YZ(i) = 0.0
c	   write(18,*) 	'Perturbed Residue', i
	   do j=1,resnum	  
c			write(18,*) i,j,delRperb_X_XYZ(i,j)
	   		sum_X_YZ(i) = sum_X_YZ(i) + delRperb_X_YZ(i,j)
	   		sum_Y_YZ(i) = sum_Y_YZ(i) + delRperb_Y_YZ(i,j)
	   		sum_Z_YZ(i) = sum_Z_YZ(i) + delRperb_Z_YZ(i,j)
	   enddo
			sum_X_YZ(i) = sum_X_YZ(i)/resnum
			sum_Y_YZ(i) = sum_Y_YZ(i)/resnum
			sum_Z_YZ(i) = sum_Z_YZ(i)/resnum
	   enddo
	   

c	   do i=1,resnum
c			write(14,*) sum_X(i),sum_Y(i),sum_Z(i)
c	   enddo


	   sumNormYZ=0.0
	   do i=1,resnum
		distS2_YZ(i) = sqrt(sum_X_YZ(i)*sum_X_YZ(i)+sum_Y_YZ(i)*sum_Y_YZ(i)+sum_Z_YZ(i)*sum_Z_YZ(i))
		sumNormYZ = sumNormYZ + distS2_YZ(i)
	   enddo
	   do i=1,resnum
		 distS2_YZ(i)= distS2_YZ(i)/sumNormYZ
	   enddo



c--- normalize matrix based on row

	   sum3=0.0
	   do i=1,resnum
	   sum2(i) = 0.0
	   do j=1,resnum
	   		sum2(i) = sum2(i) + delRperbDistMat_YZ(i,j)
	   enddo
	   sum3 = sum3 + sum2(i)
	   enddo
	   
	   do i=1,resnum
	   	do j=1,resnum
	   		delRperbDistMat_YZNorm(i,j) = delRperbDistMat_YZ(i,j) /sum2(i)
	   	enddo
	   enddo

c  column sum
	   DO i=1,RESNUM
	   sum1jyz(i) = 0.0
		DO j=1,RESNUM
			sum1jyz(i) = sum1jyz(i) + delRperbDistMat_YZ(j,i)
		ENDDO
	   ENDDO

c  ROW sum
	   DO i=1,RESNUM
	   sum1iyz(i) = 0.0
		DO j=1,RESNUM
			sum1iyz(i) = sum1iyz(i) + delRperbDistMat_YZ(i,j)
		ENDDO
	   ENDDO

	
	WRITE(69,*)'YZ direction'
	WRITE(*,511)ISAY,RMS0,PR0,ICHOSEN
	WRITE(69,511)ISAY,RMS0,PR0,ICHOSEN

	RMS0X = RMS0
	DO I=1,RESNUM
		XYZ(I) = XTEMP_YZ(I) 
		YYZ(I) = YTEMP_YZ(I) 
		ZYZ(I) = ZTEMP_YZ(I) 
	ENDDO

	WRITE(71,*)'YZ direction'
	WRITE(71,*)'STEP',ISAY
	do I=1,RESNUM	
      write(71,55)"ATOM  ",I," CA ",CNAM(I),"A",I,XYZ(I),YYZ(I),ZYZ(I)
     +	  ,RNUM(I),BETA(I)
	WRITE(71,*)I,DIFBST_YZ(I)
	enddo
	write(71,56)'END'



c********************************************************************
c********************************************************************

c 	Perturbation through  XYZ-direction
	
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO I=1,RESNUM
		Xp(I) = 0.
 		Yp(I) = 0.
		Zp(I) = 0.
		Xn(I) = 0.
		Yn(I) = 0.
		Zn(I) = 0.
	ENDDO
    
     
	
	DO 1006 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-2,3*nperturbRes
c			delForce(i) = DELF*(ran1(idum)*2-1)
			delForce(i) = 1.0
		ENDDO
		

C	Multiply force vector with inverse of Hessian = DeltaR

	DO I=1,RES3
	 delRperb_XYZ(I)= 0.0
	 DO J=1,RES3
	    delRperb_XYZ(I)= delRperb_XYZ(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
	    XI_XYZ(I) = XI(I)
        YI_XYZ(I) = YI(I)
        ZI_XYZ(I) = ZI(I)	
	ENDDO

c calculate new coordinates based on any ratio
	DO I=perbegin,perend
	   Xp(I) = XI_XYZ(I) + alpha * delRperb_XYZ(3*I-2)
	   Yp(I) = YI_XYZ(I) + alpha * delRperb_XYZ(3*I-1)
	   Zp(I) = ZI_XYZ(I) + alpha * delRperb_XYZ(3*I)
	ENDDO

	DO I=perbegin,perend
	   Xn(I) = XI_XYZ(I) - alpha * delRperb_XYZ(3*I-2)
	   Yn(I) = YI_XYZ(I) - alpha * delRperb_XYZ(3*I-1)
	   Zn(I) = ZI_XYZ(I) - alpha * delRperb_XYZ(3*I)
	ENDDO

c write into output files
c	if(nperturbRes.lt.10) then
c		write(49,411) nperturbRes
c		write(90,411) nperturbRes
c		write(81,411) nperturbRes
c	else if(nperturbRes.lt.100) then
c		write(49,412) nperturbRes
c		write(90,412) nperturbRes
c		write(81,412) nperturbRes
c	else
c		write(49,413) nperturbRes
c		write(90,413) nperturbRes
c		write(81,413) nperturbRes
c	endif
c411		FORMAT('PerturbRes',I1)
c412		FORMAT('PerturbRes',I2)
c413		FORMAT('PerturbRes',I3)

c	do I=1,RESNUM
c		write(49,*)I, Xp(I),Yp(I),Zp(I)
c		write(90,*)I, Xn(I),Yn(I),Zn(I)
c	enddo


	DO I=1,RESNUM
	    XT_XYZ(I) = 0.
        YT_XYZ(I) = 0.
        ZT_XYZ(I) = 0.	
	ENDDO

	do I=1,RESNUM	
		XT_XYZ(I) = XI_XYZ(I) +  delRperb_XYZ(3*I-2) 
		YT_XYZ(I) = YI_XYZ(I) +  delRperb_XYZ(3*I-1) 
		ZT_XYZ(I) = ZI_XYZ(I) +  delRperb_XYZ(3*I) 
	enddo

C	Superimpose the new structure to the initial,

	call impose(NR,XI_XYZ,YI_XYZ,ZI_XYZ,NR,XT_XYZ,YT_XYZ,ZT_XYZ,rmsvalue)


C	COMPUTE THE DIFFERENCES

	DO I=1,RESNUM
		DX_XYZ=XT_XYZ(I)-XI_XYZ(I)
		DY_XYZ=YT_XYZ(I)-YI_XYZ(I)
		DZ_XYZ=ZT_XYZ(I)-ZI_XYZ(I)
		DIFF_XYZ(I)=SQRT(DX_XYZ*DX_XYZ+DY_XYZ*DY_XYZ+DZ_XYZ*DZ_XYZ)
	ENDDO

	call impose(NR,XF,YF,ZF,NR,XT_XYZ,YT_XYZ,ZT_XYZ,rmsvalue)
	CALL pearsn (DIFF_XYZ,DTARGET,NR,PR)

	WRITE(72,*)'XYZ diRection'
c	WRITE(*,*)perturbRes,rmsvalue,PR
	WRITE(72,*)nperturbRes,rmsvalue,PR

	CALL INDEXX(RESNUM,PRCORR,INDX)


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
c	IF(PR.GT.PR0)THEN
c	IF(PR.GT.PR0.AND.rmsvalue.LT.RMS0)THEN
C	IF(rmsvalue.LT.RMS0)THEN
		PR0=PR
		RMS0=rmsvalue
		ICHOSEN=nperturbRes
		DO I=1,RESNUM
			XTEMP_XYZ(I) = XT_XYZ(I) 
			YTEMP_XYZ(I) = YT_XYZ(I) 
			ZTEMP_XYZ(I) = ZT_XYZ(I) 
			DIFBST_XYZ(I) = DIFF_XYZ(I)
		ENDDO
c	ENDIF

	do I=perbegin,perend
	   ratioXDiff_XYZ = alpha * (XTEMP_XYZ(I)-XI_XYZ(I))
	   ratioYDiff_XYZ = alpha * (XTEMP_XYZ(I)-YI_XYZ(I))
	   ratioZDiff_XYZ = alpha * (XTEMP_XYZ(I)-ZI_XYZ(I))
	   SQDIS2ratio_XYZ(nperturbRes,I) = sqrt ( ratioXDiff_XYZ*ratioXDiff_XYZ + ratioYDiff_XYZ * ratioYDiff_XYZ + ratioZDiff_XYZ * ratioZDiff_XYZ )
c		write(81,*)I, ratioXDiff_XYZ,ratioYDiff_XYZ,ratioZDiff_XYZ,SQDIS2ratio_XYZ(nperturbRes,I)
	enddo

c ------------------



	   sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_XYZ(i) = sqrt( delRperb_XYZ(3*I-2) *  delRperb_XYZ(3*I-2) +  delRperb_XYZ(3*I-1) * delRperb_XYZ(3*I-1) +  delRperb_XYZ(3*I) *  delRperb_XYZ(3*I) )
		   sum = sum + delRperbDist_XYZ(i) 		   
	   ENDDO
	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_XYZ(i,j) = delRperbDist_XYZ(j) 
			delRperb_X_XYZ(i,j) = delRperb_XYZ(3*J-2) 
			delRperb_Y_XYZ(i,j) = delRperb_XYZ(3*J-1) 
			delRperb_Z_XYZ(i,j) = delRperb_XYZ(3*J) 
		ENDDO
	   ENDDO

c---------------------
	  
	   
1006	enddo

c---  delRperb_X_XYZ,delRperb_Y_XYZ,delRperb_Z_XYZ wrt row



	   do i=1,resnum
	   sum_X_XYZ(i) = 0.0
	   sum_Y_XYZ(i) = 0.0
	   sum_Z_XYZ(i) = 0.0
c	   write(18,*) 	'Perturbed Residue', i
	   do j=1,resnum	  
c			write(18,*) i,j,delRperb_X_XYZ(i,j)
	   		sum_X_XYZ(i) = sum_X_XYZ(i) + delRperb_X_XYZ(i,j)
	   		sum_Y_XYZ(i) = sum_Y_XYZ(i) + delRperb_Y_XYZ(i,j)
	   		sum_Z_XYZ(i) = sum_Z_XYZ(i) + delRperb_Z_XYZ(i,j)
	   enddo
			sum_X_XYZ(i) = sum_X_XYZ(i)/resnum
			sum_Y_XYZ(i) = sum_Y_XYZ(i)/resnum
			sum_Z_XYZ(i) = sum_Z_XYZ(i)/resnum
	   enddo
	   

c	   do i=1,resnum
c			write(14,*) sum_X(i),sum_Y(i),sum_Z(i)
c	   enddo


	   sumNormXYZ=0.0
	   do i=1,resnum
		distS2_XYZ(i) = sqrt(sum_X_XYZ(i)*sum_X_XYZ(i)+sum_Y_XYZ(i)*sum_Y_XYZ(i)+sum_Z_XYZ(i)*sum_Z_XYZ(i))
		sumNormXYZ = sumNormXYZ + distS2_XYZ(i)
	   enddo
	   do i=1,resnum
		 distS2_XYZ(i)= distS2_XYZ(i)/sumNormXYZ
	   enddo


c--- normalize matrix delRperbDistMat_XYZ wrt row

	   sum3=0.0
	   do i=1,resnum
	   sum2(i) = 0.0
	   do j=1,resnum
	   		sum2(i) = sum2(i) + delRperbDistMat_XYZ(i,j)
	   enddo
	   sum3 = sum3 + sum2(i)
	   enddo
	   
	   do i=1,resnum
	   	do j=1,resnum
	   		delRperbDistMat_XYZNorm(i,j) = delRperbDistMat_XYZ(i,j) /sum2(i)
	   	enddo
	   enddo

c  column sum
	   DO i=1,RESNUM
	   sum1jxyz(i) = 0.0
		DO j=1,RESNUM
			sum1jxyz(i) = sum1jxyz(i) + delRperbDistMat_XYZ(j,i)
		ENDDO
	   ENDDO

c row sum
	   DO i=1,RESNUM
	   sum1ixyz(i) = 0.0
		DO j=1,RESNUM
			sum1ixyz(i) = sum1ixyz(i) + delRperbDistMat_XYZ(i,j)
		ENDDO
c		write(*,*) i,sum1ixyz(i)
	   ENDDO

	
	WRITE(69,*)'XYZ direction'
	WRITE(*,511)ISAY,RMS0,PR0,ICHOSEN
	WRITE(69,511)ISAY,RMS0,PR0,ICHOSEN

	RMS0X = RMS0
	DO I=1,RESNUM
		XXYZ(I) = XTEMP_XYZ(I) 
		YXYZ(I) = YTEMP_XYZ(I) 
		ZXYZ(I) = ZTEMP_XYZ(I) 
	ENDDO

	WRITE(71,*)'XYZ direction'
	WRITE(71,*)'STEP',ISAY
	do I=1,RESNUM	
c      write(71,55)"ATOM  ",I," CA ",CNAM(I),"A",I,XXYZ(I),YXYZ(I),ZXYZ(I)
c     +	  ,RNUM(I),BETA(I)
	WRITE(71,*)I,DIFBST_XYZ(I)
	enddo
	write(71,56)'END'

c-- average and max of distance from responseVector

	do i=1,resnum
		AvgdistS2(i)=(distS2_X(i)+distS2_Y(i)+distS2_Z(i)+distS2_XY(i)+distS2_XZ(i)+distS2_YZ(i)+distS2_XYZ(i))/7
		MaxdistS2(i)=max(distS2_X(i),distS2_Y(i),distS2_Z(i),distS2_XY(i),distS2_XZ(i),distS2_YZ(i),distS2_XYZ(i))
		write(15,*) distS2_XYZ(i),AvgdistS2(i)
		write(16,*) MaxdistS2(i)
	enddo

c-----  average column sum 

	sumCol = 0.0
	do i=1,resnum
		AvgCol(i)= (sum1jx(i)+sum1jy(i)+sum1jz(i)+sum1jxy(i)+sum1jyz(i)+sum1jxz(i)+sum1jxyz(i))/7
		sumCol = sumCol + AvgCol(i)
	enddo
	do i=1,resnum
		AvgnormCol(i)= AvgCol(i)/sumCol
	enddo

	do i=1,resnum
		write(88,*) AvgnormCol(i) !AvgS2i
	enddo
	
c-----  maximum column sum 
c
	sumCol1=0
	do i=1,resnum
		MaxCol(i)= Max(sum1jx(i),sum1jy(i),sum1jz(i),sum1jxy(i),sum1jyz(i),sum1jxz(i),sum1jxyz(i))
		sumCol1 = sumCol1 + MaxCol(i)
	enddo
	do i=1,resnum
		MaxnormCol(i)= MaxCol(i)/sumCol1
	enddo

	do i=1,resnum
		write(89,*) MaxnormCol(i) !MaxS2i
	enddo

c-----  average row sum 

	sumRow = 0.0
	do i=1,resnum
		AvgRow(i)= (sum1ix(i)+sum1iy(i)+sum1iz(i)+sum1ixy(i)+sum1iyz(i)+sum1ixz(i)+sum1ixyz(i))/7
		sumRow = sumRow + AvgRow(i)
	enddo
	do i=1,resnum
		AvgnormRow(i)= AvgRow(i)/sumRow
	enddo

	do i=1,resnum
		write(90,*) AvgnormRow(i) !AvgRowSum
	enddo
	
c-----  maximum row sum 
c
	sumRow1=0
	do i=1,resnum
		MaxRow(i)= Max(sum1ix(i),sum1iy(i),sum1iz(i),sum1ixy(i),sum1iyz(i),sum1ixz(i),sum1ixyz(i))
		sumRow1 = sumRow1 + MaxRow(i)
	enddo
	do i=1,resnum
		MaxnormRow(i)= MaxRow(i)/sumRow1
	enddo

	do i=1,resnum
		write(91,*) MaxnormRow(i) !MaxRowSum
	enddo



c------------------------------------------------
c apply svd to the normalized XYZ response matrix (based on magnitude)
c------------------------------------------------
c	   do i=1,resnum
c	   		do j=1,resnum
c	   delRperbDistMat_Avg(i,j) = (delRperbDistMat_XNorm(i,j)+delRperbDistMat_YNorm(i,j)+delRperbDistMat_ZNorm(i,j)+delRperbDistMat_XZNorm(i,j)+delRperbDistMat_YZNorm(i,j)+delRperbDistMat_XYZNorm(i,j)+delRperbDistMat_XYNorm(i,j))/7.0		
c			delRperbDistMat_Avg(i,j) = (delRperbDistMat_XYZNorm(i,j))
c     		enddo
c      enddo

c		do i=1,resnum
c				write(82,128) (delRperbDistMat_Avg(i,j),j=1,resnum)
c		enddo
c128	format(500f8.5)

c	    CALL SVDCMP(delRperbDistMat_Avg,resnum,resnum,resnum,resnum,Wr,Vr)
c	    CALL INDEXX(resnum,Wr,INDX)
c	    do i=1,resnum
c	    	write(*,*) i,INDX(i),Wr(INDX(i))
c	    	write(35,*) i,INDX(resnum),Wr(INDX(resnum)),(Vr(i,INDX(k)),k=resnum,resnum,1),(delRperbDistMat_Avg(i,INDX(k)),k=resnum,resnum,1)
c	    enddo


c	do i=1,resnum
c		u12(i) =Vr(i,INDX(resnum))*Vr(i,INDX(resnum))
c		u22(i) =Vr(i,INDX(resnum-1))*Vr(i,INDX(resnum-1))
c		write(89,*) u12(i),u22(i) !pca1
c	enddo
	



	WRITE(*,*) 'Program finished successfully! Please press enter.'

666	STOP
	END


C	THE SUPERPOSITION SUBROUTINE IS ADAPTED FROM THE TINKER PROGRAM
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine impose  --  superimpose two coordinate sets  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "impose" performs the least squares best superposition
c     of two atomic coordinate sets via a quaternion method;
c     upon return, the first coordinate set is unchanged while
c     the second set is translated and rotated to give best fit;
c     the final root mean square fit is returned in "rmsvalue"
c
c
	subroutine impose (n1,x1,y1,z1,n2,x2,y2,z2,rmsvalue)
	IMPLICIT REAL*8 (A-H,O-Z)
	integer n1
      REAL*8 x1(n1),y1(n1),z1(n1)
      REAL*8 x2(n2),y2(n2),z2(n2)
      REAL*8 xmid,ymid,zmid
      REAL*8 rmsvalue,rmsfit
      REAL*8 wfit(n1)
	dimension ifit(2,n1)
	logical verbose
	integer nfit
c
c
c     superimpose the full structures if not specified%
c	superposition will best-fit on atoms for wfit=1.
c


	nfit = 0  
	verbose=.FALSE.
      if (nfit .eq. 0) then
         nfit = min(n1,n2)
         do i = 1, nfit
            ifit(1,i) = i
            ifit(2,i) = i
            wfit(i) = 1.0d0
c            wfit(i) = 0.0d0
c		  if((i.ge.83.and.i.le.87).or.(i.ge.102.and.i.le.225).or.
c     +	  (i.ge.277.and.i.le.307)) wfit(i) = 1.0d0
         end do
      end if
c
c     find the rms fit of input coordinates
c



      if (verbose) then
         rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2,n1,n2,ifit,wfit,nfit)
         write (*,20)  rmsvalue
   20    format (/,' IMPOSE  --  Input Coordinates',12x,f12.6)
      end if
c

c     superimpose the centroids of active atom pairs
c


	  call centr (n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid,ifit,wfit,nfit)
	  
	  
	  
      if (verbose) then
         rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2,n1,n2,ifit,wfit,nfit)
         write (*,30)  rmsvalue
   30    format (' IMPOSE  --  After Translation',12x,f12.6)
      end if
c
c     use a quaternion method to achieve the superposition
c
      call quatfit (n1,x1,y1,z1,n2,x2,y2,z2,ifit,wfit,nfit)

      rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2,n1,n2,ifit,wfit,nfit)
      if (verbose) then
         write (*,40)  rmsvalue
   40    format (' IMPOSE  --  After Rotation',15x,f12.6)
      end if
c
c     translate both coordinate sets so as to return
c     the first set to its original position
c
      do i = 1, n1
         x1(i) = x1(i) + xmid
         y1(i) = y1(i) + ymid
         z1(i) = z1(i) + zmid
      end do
      do i = 1, n2
         x2(i) = x2(i) + xmid
         y2(i) = y2(i) + ymid
         z2(i) = z2(i) + zmid
      end do

	  nfit=0
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine center  --  superimpose structure centroids  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "centr" moves the weighted centroid of each coordinate
c     set to the origin during least squares superposition
c
c
      subroutine centr(n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid,ifit,wfit
     +,nfit)
      IMPLICIT REAL*8 (A-H,O-Z)
	integer i,k,n1,n2
      REAL*8 xmid,ymid,zmid,weight,norm
      REAL*8 x1(n1),y1(n1),z1(n1)
      REAL*8 x2(n2),y2(n2),z2(n2)
      REAL*8 wfit(n1)
	dimension ifit(2,n1)
c
c
c     find the weighted centroid of the second
c     structure and translate it to the origin
c
	    
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(2,i)
         weight = wfit(i)
         xmid = xmid + x2(k)*weight
         ymid = ymid + y2(k)*weight
         zmid = zmid + z2(k)*weight
         norm = norm + weight
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n2
         x2(i) = x2(i) - xmid
         y2(i) = y2(i) - ymid
         z2(i) = z2(i) - zmid
      end do
c
c     now repeat for the first structure, note
c     that this centroid position gets returned
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(1,i)
         weight = wfit(i)
         xmid = xmid + x1(k)*weight
         ymid = ymid + y1(k)*weight
         zmid = zmid + z1(k)*weight
         norm = norm + weight
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n1
         x1(i) = x1(i) - xmid
         y1(i) = y1(i) - ymid
         z1(i) = z1(i) - zmid
      end do
	  
	  
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine quatfit  --  quaternion superposition of coords  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "quatfit" uses a quaternion-based method to achieve the best
c     fit superposition of two sets of coordinates
c
c     literature reference:
c
c     S. J. Kearsley, "An Algorithm for the Simultaneous Superposition
c     of a Structural Series", Journal of Computational Chemistry,
c     11, 1187-1192 (1990)
c
c     adapted from an original program written by David J. Heisterberg,
c     Ohio Supercomputer Center, Columbus, OH
c
c
      subroutine quatfit (n1,x1,y1,z1,n2,x2,y2,z2,ifit,wfit,nfit)
      IMPLICIT REAL*8 (A-H,O-Z)
	integer i,i1,i2,n1,n2
      REAL*8 weight,xrot,yrot,zrot
      REAL*8 xxyx,xxyy,xxyz,xyyx,xyyy
      REAL*8 xyyz,xzyx,xzyy,xzyz
      REAL*8 rot(3,3),temp1(4),temp2(4)
      REAL*8 q(4),d(4),c(4,4),v(4,4)
      REAL*8 x1(n1),y1(n1),z1(n1)
      REAL*8 x2(n2),y2(n2),z2(n2)
      REAL*8 wfit(n1)
	dimension ifit(2,n1)

c
c
c     build the upper triangle of the quadratic form matrix
c


      xxyx = 0.0d0
      xxyy = 0.0d0
      xxyz = 0.0d0
      xyyx = 0.0d0
      xyyy = 0.0d0
      xyyz = 0.0d0
      xzyx = 0.0d0
      xzyy = 0.0d0
      xzyz = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weight = wfit(i)
         xxyx = xxyx + weight*x1(i1)*x2(i2)
         xxyy = xxyy + weight*y1(i1)*x2(i2)
         xxyz = xxyz + weight*z1(i1)*x2(i2)
         xyyx = xyyx + weight*x1(i1)*y2(i2)
         xyyy = xyyy + weight*y1(i1)*y2(i2)
         xyyz = xyyz + weight*z1(i1)*y2(i2)
         xzyx = xzyx + weight*x1(i1)*z2(i2)
         xzyy = xzyy + weight*y1(i1)*z2(i2)
         xzyz = xzyz + weight*z1(i1)*z2(i2)
      end do
      c(1,1) = xxyx + xyyy + xzyz
      c(1,2) = xzyy - xyyz
      c(2,2) = xxyx - xyyy - xzyz
      c(1,3) = xxyz - xzyx
      c(2,3) = xxyy + xyyx
      c(3,3) = xyyy - xzyz - xxyx
      c(1,4) = xyyx - xxyy
      c(2,4) = xzyx + xxyz
      c(3,4) = xyyz + xzyy
      c(4,4) = xzyz - xxyx - xyyy
c
c     diagonalize the quadratic form matrix
c
      call jacobi (4,4,c,d,v,temp1,temp2)
c
c     extract the desired quaternion
c
      q(1) = v(1,4)
      q(2) = v(2,4)
      q(3) = v(3,4)
      q(4) = v(4,4)
c
c     assemble rotation matrix that superimposes the molecules
c
      rot(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
      rot(2,1) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
      rot(3,1) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
      rot(1,2) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
      rot(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
      rot(3,2) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
      rot(1,3) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
      rot(2,3) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
      rot(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
c
c     rotate second molecule to best fit with first molecule
c
      do i = 1, n2
         xrot = x2(i)*rot(1,1) + y2(i)*rot(1,2) + z2(i)*rot(1,3)
         yrot = x2(i)*rot(2,1) + y2(i)*rot(2,2) + z2(i)*rot(2,3)
         zrot = x2(i)*rot(3,1) + y2(i)*rot(3,2) + z2(i)*rot(3,3)
         x2(i) = xrot
         y2(i) = yrot
         z2(i) = zrot
      end do
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine jacobi  --  jacobi matrix diagonalization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "jacobi" performs a matrix diagonalization of a real
c     symmetric matrix by the method of Jacobi rotations
c
c     n    logical dimension of the matrix to be diagonalized
c     np   physical dimension of the matrix storage area
c     a    input with the matrix to be diagonalized; only
c             the upper triangle and diagonal are required
c     d    returned with the eigenvalues in ascending order
c     v    returned with the eigenvectors of the matrix
c     b    temporary work vector
c     z    temporary work vector
c
c

      subroutine jacobi (n,np,a,d,v,b,z)
	IMPLICIT REAL*8 (A-H,O-Z)
      integer i,j,k,ip,iq,n,np,nrot,maxrot
      REAL*8 sm,tresh,s,c,t,theta,tau,h,g,p
      REAL*8 a(np,np),d(np),v(np,np),b(np),z(np)
c
c
c     setup and initialization
c
      maxrot = 100
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
c
c     perform the jacobi rotations
c
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip))
     &                    .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
c
c     print warning if not converged
c
   10 continue
      if (nrot .eq. maxrot) then
         write (*,20)
   20    format (/,' JACOBI  --  Matrix Diagonalization not Converged')
      end if
c
c     sort the eigenvalues and vectors
c
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end

      function rmsfit (x1,y1,z1,x2,y2,z2,n1,n2,ifit,wfit,nfit)
      IMPLICIT REAL*8 (A-H,O-Z)
	integer i,i1,i2,n1
      REAL*8 rmsfit,rmsterm
      REAL*8 xr,yr,zr,dist2
      REAL*8 weight,norm
      REAL*8 x1(n1),y1(n1),z1(n1)
      REAL*8 x2(n2),y2(n2),z2(n2)
      REAL*8 wfit(n1)
	dimension ifit(2,n1)
c
c
c     compute the rms fit over superimposed atom pairs
c

      rmsfit = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weight = wfit(i)
         xr = x1(i1) - x2(i2)
         yr = y1(i1) - y2(i2)
         zr = z1(i1) - z2(i2)
         dist2 = xr**2 + yr**2 + zr**2
         norm = norm + weight
         rmsterm = dist2 * weight
         rmsfit = rmsfit + rmsterm
      end do
      rmsfit = sqrt(rmsfit/norm)
      return
      end


C*****************************************************************
C	ALL OF THE FOLLOWING FUNCTIONS ARE TAKEN FROM 
C     PRESS, W.H. ET AL. "NUMERICAL RECIPES IN FORTRAN 77", 
C	CAMBRIDGE UNIVERSITY PRESS, 2001
C*****************************************************************
C															   
C	SINGULAR VALUE DECOMPOSITION							   
C															   
C*****************************************************************

      SUBROUTINE SVDCMP(a,m,n,mp,np,w,v)
	IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=70000)
C     USES pythag
      INTEGER i,its,j,jj,k,l,nm
      REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0
33      continue
        if(g.ne.0.0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.D0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END

      FUNCTION pythag(a,b)
      IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 a,b,pythag
      REAL*8 absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END
     

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      IF(N.EQ.1)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END




      SUBROUTINE pearsn(x,y,n,r)
c	Pearson correlation between two data sets
	IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER n
      REAL*8 prob,r,z,x(n),y(n),TINY
      PARAMETER (TINY=1.e-20)
CU    USES betai
      INTEGER j
      REAL*8 ax,ay,df,sxx,sxy,syy,t,xt,yt,betai
      ax=0.
      ay=0.
      do 11 j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
11    continue
      ax=ax/n
      ay=ay/n
      sxx=0.
      syy=0.
      sxy=0.
      do 12 j=1,n
        xt=x(j)-ax
        yt=y(j)-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
12    continue
      r=sxy/(sqrt(sxx*syy)+TINY)
      return
      END


      FUNCTION ran1(idum)
	IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
