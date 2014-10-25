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
C############################################################
	PARAMETER (NR=836,GAMMA=10.0,nbG=1)
C############################################################
	REAL*8 X(NR),Y(NR),Z(NR),XF(NR),YF(NR),ZF(NR)

	CHARACTER CNAM(NR)*3,CNAMF(NR)*3
	REAL*8 BETA(NR),HBETA(NR),conNum(NR),BETAF(NR)	
	REAL*8 BKBT(NR*3,NR*3)
	REAL*8 HESS(NR*3,NR*3),INVBKBT(NR*3,NR*3)
	REAL*8 W(NR*3),V(NR*3,NR*3)
	DIMENSION INDX(NR*3)
	REAL*8 INVCONT(NR,NR)
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
	REAL*8 delForce(NR*3),delRperb(NR*3),delR(NR*3,NR*3),delRf(NR*3),delrs(NR*3)
	INTEGER nperturbRes,iperturbRes,fperturbRes
	REAL*8 XT(NR),YT(NR),ZT(NR),DIFF(NR),DTARGET(NR)
	REAL*8 XI(NR),YI(NR),ZI(NR)
	REAL*8 XTEMP(NR),YTEMP(NR),ZTEMP(NR),DIFBST(NR),rmsvalue,PR                
	REAL*8 delRperb_X(NR*3),delRperb_Y(NR*3),delRperb_Z(NR*3)
	REAL*8 delRperb_XY(NR*3),delRperb_YZ(NR*3),delRperb_XZ(NR*3)
	REAL*8 delRperb_XYZ(NR*3)

	REAL*8 PRCORR(NR),INVBKBTRes42(3,3),INVBKBTRes36(3,3),v36(3,3),v42(3,3),INVBKBTRes(3,3),vRes(3,3),wRes(3),traceInvBKBT(NR)

	REAL*8 Xp(NR),Yp(NR),Zp(NR),Xn(NR),Yn(NR),Zn(NR)
	REAL*8 SQDIS2ratio_X(NR,NR), SQDIS2ratio_Y(NR,NR), SQDIS2ratio_Z(NR,NR)
	REAL*8 SQDIS2ratio_XY(NR,NR), SQDIS2ratio_XZ(NR,NR), SQDIS2ratio_YZ(NR,NR)
	REAL*8 SQDIS2ratio_XYZ(NR,NR), AVG(NR,NR), sumDIS2RATIO(NR)


C	DUMMIES

	INTEGER DINT,DINT2,ICA,IDUM
	REAL*8 DIFXX,DIFYY,DIFFZZ,DIST,DSUM,DSUM1,DSUM2,TSvibSum
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


	REAL*8 delRperbDist_X(NR),delRperbDistMat_X(NR,NR)
	REAL*8 delRperbDist_Y(NR),delRperbDistMat_Y(NR,NR)
	REAL*8 delRperbDist_Z(NR),delRperbDistMat_Z(NR,NR)
	REAL*8 delRperbDist_XY(NR),delRperbDistMat_XY(NR,NR)
	REAL*8 delRperbDist_XZ(NR),delRperbDistMat_XZ(NR,NR)
	REAL*8 delRperbDist_YZ(NR),delRperbDistMat_YZ(NR,NR)
	REAL*8 delRperbDist_XYZ(NR),delRperbDistMat_XYZ(NR,NR)
	
	REAL*8 bindSumRow_X(NR),allSumColumn_X(NR)
	REAL*8 bindSumRow_Y(NR),allSumColumn_Y(NR)
	REAL*8 bindSumRow_Z(NR),allSumColumn_Z(NR)
	REAL*8 bindSumRow_XY(NR),allSumColumn_XY(NR)
	REAL*8 bindSumRow_XZ(NR),allSumColumn_XZ(NR)
	REAL*8 bindSumRow_YZ(NR),allSumColumn_YZ(NR)
	REAL*8 bindSumRow_XYZ(NR),allSumColumn_XYZ(NR)


C******************************************************************
C	PARAMETERS
C******************************************************************

	RESNUM=NR
	RES3=NR*3
        RESNUM2=RESNUM/4
        WRITE(*,*) 'NUMBER OF RESNUM2',RESNUM2
        
C	radius threshold to decide if two residues/nucleotides are 
C	connected
	OPEN(15,FILE='input.parms')
C############################################################
c	READ(15,*)iperturbRes,fperturbRes,noPeptideSize,iBP1,iBP2,iBP3,iBP4,iBP5,iBP6
	READ(15,*)iperturbRes,fperturbRes,noPeptideSize
	EIGENCUT=1E-6

C############################################################


C******************************************************************
C	FILES
C******************************************************************

C	THESE ARE THE INPUT FILES
C	files can be obtained from Brookhaven Protein Data Bank 
C	http://www.rcsb.org/pdb
C	the 1st is the structure to be perturbed, the 2nd is the target 
C############################################################
	OPEN(50,FILE='4DXMmod_atmorder.pdb')
	OPEN(51,FILE='4DXMmod_atmorder.pdb')
C############################################################

c *******output files


	OPEN(60,FILE='init_centers.pdb')
	OPEN(59,FILE='target_centers.pdb')
	open(61,file='4DXM.delRinitial')
	OPEN(66,FILE='eigenvalues.txt')
	open(69,file='4DXM_RESULTS.txt')
	OPEN(71,FILE='4DXM_diffbest.txt')
	OPEN(72,FILE='4DXM_dev_corr.txt')


	OPEN (95,file='20sloweigenx.txt')
	OPEN (96,file='20sloweigeny.txt')
	OPEN (97,file='20sloweigenz.txt')
	OPEN (98,file='20sloweigenr.txt')

	
c ********************
c X -direction 
c ********************
	open(27,file='4DXM_X_BindSumRows.dat')
	open(28,file='4DXM_X_AllSumColumn.dat')

c ********************
c Y -direction 
c ********************
	open(31,file='4DXM_Y_BindSumRows.dat')
	open(32,file='4DXM_Y_AllSumColumn.dat')

c ********************
c Z -direction 
c ********************
	open(35,file='4DXM_Z_BindSumRows.dat')
	open(36,file='4DXM_Z_AllSumColumn.dat')


c ********************
c  XY -direction 
c ********************
	open(39,file='4DXM_XY_BindSumRows.dat')
	open(40,file='4DXM_XY_AllSumColumn.dat')

c ********************
c  XZ -direction 
c ********************
	open(43,file='4DXM_XZ_BindSumRows.dat')
	open(44,file='4DXM_XZ_AllSumColumn.dat')

c ********************
c  YZ -direction 
c ********************
	open(47,file='4DXM_YZ_BindSumRows.dat')
	open(48,file='4DXM_YZ_AllSumColumn.dat')

c ********************
c  XYZ -direction 
c ********************
	open(54,file='4DXM_XYZ_BindSumRows.dat')
	open(55,file='4DXM_XYZ_AllSumColumn.dat')

	open(11,file='4DXM.BindRatio_all_directions.txt')
	open(12,file='4DXM_MaxAR.txt')



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
C	one has the option of pertubing in a non-random manner; see DELF loop below
C	PRTARGET IS THE TARGETED PEARSON CORRELATION
C	RMSTARGET IS THE TARGETED RMSD BETWEEN THE PERTURBED STRUCTURE AND THE EXPERIMENTAL.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	ITERF=0
	DELF=0.1
	PRTARGET=0.999
	RMSTARGET=0.1
c	subroutine ran1 uses needs a negative integer seed for initialization:
c	idum=-67126379

500	CONTINUE

	ISAY=ISAY+1

	write(72,*)"STEP ",ISAY
	write(66,*)"STEP ",ISAY

C	PR0 will store the highest Pearson correlation in each iteration
C	RMS0 will store the lowest RMSD in each iteration

	PR0=0.
	RMS0=1000000.

C******************************************************************
C	construct directional cosine matrix (B)
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
	
	 WRITE(*,*) 'NO. OF MEMBERS:',totConMax


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
C	MODESHAPES - Twenty Slowest Modes
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
C	IS SELECTED AS THE NEXT CANDIDATE TO BE PERTURBED
C******************************************************************






c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

c                            X-direction
	
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	
	DO 1000 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
		DO i=3*nperturbRes-2,3*nperturbRes-2
c			delForce(i) = DELF*(ran1(idum)*2-1)
c			delForceRand(i) = delForce(i)
			delForce(i) = 1.0
c		write(*,*) nperturbRes,i,delForce(i)
		ENDDO
		

C******************************************************************
C	Multiply force vector with inverse of BBT = DeltaR
C******************************************************************
c	 DO J=1,RES3
c		write(99,*) j,delForce(j)
c	 ENDDO
	 
	DO I=1,RES3
	 delRperb_X(I)= 0.0
	 DO J=1,RES3
	    delRperb_X(I)= delRperb_X(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO

     

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

C******************************************************************
C	Superimpose the new structure to the initial,
C******************************************************************

	DO I=1,RESNUM
	    XI_X(I) = X(I)
        YI_X(I) = Y(I)
        ZI_X(I) = Z(I)	
	ENDDO


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


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
	IF(PR.GT.PR0)THEN
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
	ENDIF
	
       sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_X(i) = sqrt( delRperb_X(3*I-2) *  delRperb_X(3*I-2) +  delRperb_X(3*I-1) * delRperb_X(3*I-1) +  delRperb_X(3*I) *  delRperb_X(3*I) )
c		   write(*,*) nperturbRes,i,delRperbDist(i)
		   sum = sum + delRperbDist_X(i) 		   
	   ENDDO
c	   write(*,*) sum

	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_X(i,j) = delRperbDist_X(j) / sum
		ENDDO
	   ENDDO



1000	enddo
	
c---- diagonal terms are equal to zero
	   DO I=1,RESNUM
		delRperbDistMat_X(i,i) = 0.0
	   ENDDO

c---- close neighbor interactions (i,i+1),(i,i+2),(i+3) and reverse are equal to zero
	   DO I=1,RESNUM
		delRperbDistMat_X(i,i+1) = 0.0
		delRperbDistMat_X(i,i+2) = 0.0
		delRperbDistMat_X(i,i+3) = 0.0
		delRperbDistMat_X(i+1,i) = 0.0
		delRperbDistMat_X(i+2,i) = 0.0
		delRperbDistMat_X(i+3,i) = 0.0
	   ENDDO



c this part is redundant for allosteric response
	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_X(i,j)
		ENDDO
	   ENDDO


c-----total response of residue j upon perturbations on all residues in X-direction
c--- this is the denominator of the Equation 10 in Equation 10 in Gerek and Ozkan, Plos Comp Biol,2010

	   DO i=1,RESNUM2
	   sum1 = 0.0
		DO j=1,RESNUM2
			sum1 = sum1 + delRperbDistMat_X(j,i)
		ENDDO
		write(28,*) i, sum1	
		allSumColumn_X(i)=sum1/(RESNUM2-1)
	   ENDDO


C####### we are building a matrix based on binding sites


	   DO J=1,RESNUM
	     DO I=1,52
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=54,61
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=63,66
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=68,96
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=99,108
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=110,146
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=148,149
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=151,186
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=188,208
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=210,261
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=263,270
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=272,275
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=277,305
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=308,317
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=319,355
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=357,358
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=360,395
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=397,417
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=419,470
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=472,479
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=481,484
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=486,514
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=517,526
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=528,564
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=566,567
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=569,604
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=606,626
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=628,679
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=681,688
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=690,693
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=695,723
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=726,735
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=737,773
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=775,776
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=778,813
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	     DO I=815,835
			delRperbDistMat_X(i,j) = 0.0
	  	 ENDDO
	   ENDDO	


C############################################################

c- this part is the numerator part of Equation 10 in Gerek and Ozkan, Plos Comp Biol, 20011 paper
	   DO i=1,RESNUM2
	   sum = 0.0
		DO j=1,RESNUM2
			sum = sum + delRperbDistMat_X(j,i)
		ENDDO
		write(27,*) i, sum	
		bindSumRow_X(i) = sum/noPeptideSize*4
	   ENDDO

c---- this is the redundant part
	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_X(j,i)
		ENDDO
	   ENDDO


c********************************************************************
c********************************************************************

c                            Y-direction
	
c********************************************************************
c********************************************************************



	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO 1001 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO

C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
		DO i=3*nperturbRes-1,3*nperturbRes-1
c			delForce(i) = DELF*(ran1(idum)*2-1)
c			delForceRand(i) = delForce(i)
			delForce(i) = 1.0
c		write(*,*) nperturbRes,i,delForce(i)
		ENDDO

C******************************************************************
C	Multiply force vector with inverse of BBT = DeltaR
C******************************************************************
c	 DO J=1,RES3
c		write(99,*) j,delForce(j)
c	 ENDDO
	 
	DO I=1,RES3
	 delRperb_Y(I)= 0.0
	 DO J=1,RES3
	    delRperb_Y(I)= delRperb_Y(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO

	DO I=1,RESNUM
	    XI_Y(I) = X(I)
        YI_Y(I) = Y(I)
        ZI_Y(I) = Z(I)	
	ENDDO
   

 

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

C******************************************************************
C	Superimpose the new structure to the initial,
C******************************************************************

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


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
	IF(PR.GT.PR0)THEN
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
	ENDIF

       sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_Y(i) = sqrt( delRperb_Y(3*I-2) *  delRperb_Y(3*I-2) +  delRperb_Y(3*I-1) * delRperb_Y(3*I-1) +  delRperb_Y(3*I) *  delRperb_Y(3*I) )
c		   write(*,*) nperturbRes,i,delRperbDist(i)
		   sum = sum + delRperbDist_Y(i) 		   
	   ENDDO
c	   write(*,*) sum

	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_Y(i,j) = delRperbDist_Y(j) / sum
		ENDDO
	   ENDDO


1001	enddo

	   DO I=1,RESNUM
		delRperbDistMat_Y(i,i) = 0.0
	   ENDDO


	   DO I=1,RESNUM
		delRperbDistMat_Y(i,i+1) = 0.0
		delRperbDistMat_Y(i,i+2) = 0.0
		delRperbDistMat_Y(i,i+3) = 0.0
		delRperbDistMat_Y(i+1,i) = 0.0
		delRperbDistMat_Y(i+2,i) = 0.0
		delRperbDistMat_Y(i+3,i) = 0.0
	   ENDDO


	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_Y(i,j)
		ENDDO
	   ENDDO

	   DO i=1,RESNUM2
	   sum1 = 0.0
		DO j=1,RESNUM2
			sum1 = sum1 + delRperbDistMat_Y(j,i)
		ENDDO
		write(32,*) i, sum1	
		allSumColumn_Y(i)=sum1/(RESNUM2-1)
	   ENDDO



C############################################################

	   DO J=1,RESNUM
	     DO I=1,52
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=54,61
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=63,66
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=68,96
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=99,108
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=110,146
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=148,149
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=151,186
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=188,208
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=210,261
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=263,270
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=272,275
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=277,305
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=308,317
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=319,355
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=357,358
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=360,395
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=397,417
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=419,470
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=472,479
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=481,484
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=486,514
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=517,526
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=528,564
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=566,567
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=569,604
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=606,626
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=628,679
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=681,688
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=690,693
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=695,723
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=726,735
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=737,773
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=775,776
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=778,813
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	     DO I=815,835
			delRperbDistMat_Y(i,j) = 0.0
	  	 ENDDO
	   ENDDO	

	   

C############################################################

	   DO i=1,RESNUM2
	   sum = 0.0
		DO j=1,RESNUM2
			sum = sum + delRperbDistMat_Y(j,i)
		ENDDO
		write(31,*) i, sum	
		bindSumRow_Y(i) = sum/noPeptideSize*4
	   ENDDO

	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_Y(j,i)
		ENDDO
	   ENDDO



c********************************************************************
c********************************************************************

c                            Z-direction
	
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

	DO 1002 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes,3*nperturbRes
C			delForce(i) = DELF*(ran1(idum)*2-1)
C			delForceRand(i) = delForce(i)
			delForce(i) = 1.0
c		write(*,*) nperturbRes,i,delForce(i)
		ENDDO
		

C******************************************************************
C	Multiply force vector with inverse of BBT = DeltaR
C******************************************************************
c	 DO J=1,RES3
c		write(99,*) j,delForce(j)
c	 ENDDO
	 
	DO I=1,RES3
	 delRperb_Z(I)= 0.0
	 DO J=1,RES3
	    delRperb_Z(I)= delRperb_Z(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO

	DO I=1,RESNUM
		XI_Z(I) = X(I)
        YI_Z(I) = Y(I)
        ZI_Z(I) = Z(I)	
	ENDDO
 
 
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

C******************************************************************
C	Superimpose the new structure to the initial,
C******************************************************************

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



C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
	IF(PR.GT.PR0)THEN
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
	ENDIF

       sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_Z(i) = sqrt( delRperb_Z(3*I-2) *  delRperb_Z(3*I-2) +  delRperb_Z(3*I-1) * delRperb_Z(3*I-1) +  delRperb_Z(3*I) *  delRperb_Z(3*I) )
c		   write(*,*) nperturbRes,i,delRperbDist(i)
		   sum = sum + delRperbDist_Z(i) 		   
	   ENDDO
c	   write(*,*) sum

	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_Z(i,j) = delRperbDist_Z(j) / sum
		ENDDO
	   ENDDO

1002	enddo
	
	   DO I=1,RESNUM
		delRperbDistMat_Z(i,i) = 0.0
	   ENDDO


	   DO I=1,RESNUM
		delRperbDistMat_Z(i,i+1) = 0.0
		delRperbDistMat_Z(i,i+2) = 0.0
		delRperbDistMat_Z(i,i+3) = 0.0
		delRperbDistMat_Z(i+1,i) = 0.0
		delRperbDistMat_Z(i+2,i) = 0.0
		delRperbDistMat_Z(i+3,i) = 0.0
	   ENDDO


	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_Z(i,j)
		ENDDO
	   ENDDO

	   DO i=1,RESNUM2
	   sum1 = 0.0
		DO j=1,RESNUM2
			sum1 = sum1 + delRperbDistMat_Z(j,i)
		ENDDO
		write(36,*) i, sum1	
		allSumColumn_Z(i)=sum1/(RESNUM2-1)
	   ENDDO


C############################################################

	   DO J=1,RESNUM
	     DO I=1,52
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=54,61
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=63,66
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=68,96
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=99,108
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=110,146
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=148,149
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=151,186
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=188,208
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=210,261
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=263,270
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=272,275
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=277,305
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=308,317
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=319,355
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=357,358
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=360,395
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=397,417
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=419,470
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=472,479
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=481,484
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=486,514
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=517,526
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=528,564
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=566,567
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=569,604
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=606,626
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=628,679
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=681,688
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=690,693
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=695,723
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=726,735
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=737,773
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=775,776
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=778,813
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	     DO I=815,835
			delRperbDistMat_Z(i,j) = 0.0
	  	 ENDDO
	   ENDDO	



C############################################################

	   DO i=1,RESNUM2
	   sum = 0.0
		DO j=1,RESNUM2
			sum = sum + delRperbDistMat_Z(j,i)
		ENDDO
		write(35,*) i, sum	
		bindSumRow_Z(i) = sum/noPeptideSize*4
	   ENDDO

	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_Z(j,i)
		ENDDO
	   ENDDO
	
c********************************************************************
c********************************************************************

c                            XY-direction
	
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

   
	DO 1003 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-2,3*nperturbRes-1
C			delForce(i) = DELF*(ran1(idum)*2-1)
C			delForceRand(i) = delForce(i)
			delForce(i) = 1.0
c		write(*,*) nperturbRes,i,delForce(i)
		ENDDO
		

C******************************************************************
C	Multiply force vector with inverse of BBT = DeltaR
C******************************************************************
c	 DO J=1,RES3
c		write(99,*) j,delForce(j)
c	 ENDDO
	 
	DO I=1,RES3
	 delRperb_XY(I)= 0.0
	 DO J=1,RES3
	    delRperb_XY(I)= delRperb_XY(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
	    XI_XY(I) = X(I)
        YI_XY(I) = Y(I)
        ZI_XY(I) = Z(I)	
	ENDDO



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

C******************************************************************
C	Superimpose the new structure to the initial,
C******************************************************************

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


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
	IF(PR.GT.PR0)THEN
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
	ENDIF

       sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_XY(i) = sqrt( delRperb_XY(3*I-2) *  delRperb_XY(3*I-2) +  delRperb_XY(3*I-1) * delRperb_XY(3*I-1) +  delRperb_XY(3*I) *  delRperb_XY(3*I) )
c		   write(*,*) nperturbRes,i,delRperbDist(i)
		   sum = sum + delRperbDist_XY(i) 		   
	   ENDDO
c	   write(*,*) sum

	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_XY(i,j) = delRperbDist_XY(j) / sum
		ENDDO
	   ENDDO

1003	enddo
	
	   DO I=1,RESNUM
		delRperbDistMat_XY(i,i) = 0.0
	   ENDDO


	   DO I=1,RESNUM
		delRperbDistMat_XY(i,i+1) = 0.0
		delRperbDistMat_XY(i,i+2) = 0.0
		delRperbDistMat_XY(i,i+3) = 0.0
		delRperbDistMat_XY(i+1,i) = 0.0
		delRperbDistMat_XY(i+2,i) = 0.0
		delRperbDistMat_XY(i+3,i) = 0.0
	   ENDDO



	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_XY(i,j)
		ENDDO
	   ENDDO

	   DO i=1,RESNUM2
	   sum1 = 0.0
		DO j=1,RESNUM2
			sum1 = sum1 + delRperbDistMat_XY(j,i)
		ENDDO
		write(40,*) i, sum1	
		allSumColumn_XY(i)=sum1/(RESNUM2-1)
	   ENDDO


C############################################################

	   DO J=1,RESNUM
	     DO I=1,52
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=54,61
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=63,66
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=68,96
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=99,108
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=110,146
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=148,149
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=151,186
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=188,208
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=210,261
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=263,270
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=272,275
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=277,305
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=308,317
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=319,355
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=357,358
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=360,395
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=397,417
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=419,470
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=472,479
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=481,484
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=486,514
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=517,526
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=528,564
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=566,567
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=569,604
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=606,626
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=628,679
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=681,688
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=690,693
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=695,723
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=726,735
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=737,773
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=775,776
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=778,813
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	     DO I=815,835
			delRperbDistMat_XY(i,j) = 0.0
	  	 ENDDO
	   ENDDO	
	   
C############################################################

	   DO i=1,RESNUM2
	   sum = 0.0
		DO j=1,RESNUM2
			sum = sum + delRperbDistMat_XY(j,i)
		ENDDO
		write(39,*) i, sum	
		bindSumRow_XY(i) = sum/noPeptideSize*4
	   ENDDO

	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_XY(j,i)
		ENDDO
	   ENDDO

c********************************************************************
c********************************************************************

c                            XZ-direction
	
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO
    
	DO 1004 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-2,3*nperturbRes,2
C			delForce(i) = DELF*(ran1(idum)*2-1)
C			delForceRand(i) = delForce(i)
			delForce(i) = 1.0
c		write(*,*) nperturbRes,i,delForce(i)
		ENDDO
		

C******************************************************************
C	Multiply force vector with inverse of BBT = DeltaR
C******************************************************************
c	 DO J=1,RES3
c		write(99,*) j,delForce(j)
c	 ENDDO
	 
	DO I=1,RES3
	 delRperb_XZ(I)= 0.0
	 DO J=1,RES3
	    delRperb_XZ(I)= delRperb_XZ(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
		XI_XZ(I) = X(I)
        YI_XZ(I) = Y(I)
        ZI_XZ(I) = Z(I)	
	ENDDO


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

C******************************************************************
C	Superimpose the new structure to the initial,
C******************************************************************

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
	IF(PR.GT.PR0)THEN
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
	ENDIF

       sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_XZ(i) = sqrt( delRperb_XZ(3*I-2) *  delRperb_XZ(3*I-2) +  delRperb_XZ(3*I-1) * delRperb_XZ(3*I-1) +  delRperb_XZ(3*I) *  delRperb_XZ(3*I) )
c		   write(*,*) nperturbRes,i,delRperbDist(i)
		   sum = sum + delRperbDist_XZ(i) 		   
	   ENDDO
c	   write(*,*) sum

	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_XZ(i,j) = delRperbDist_XZ(j) / sum
		ENDDO
	   ENDDO

1004	enddo
	
	   DO I=1,RESNUM
		delRperbDistMat_XZ(i,i) = 0.0
	   ENDDO


	   DO I=1,RESNUM
		delRperbDistMat_XZ(i,i+1) = 0.0
		delRperbDistMat_XZ(i,i+2) = 0.0
		delRperbDistMat_XZ(i,i+3) = 0.0
		delRperbDistMat_XZ(i+1,i) = 0.0
		delRperbDistMat_XZ(i+2,i) = 0.0
		delRperbDistMat_XZ(i+3,i) = 0.0
	   ENDDO



	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_XZ(i,j)
		ENDDO
	   ENDDO

	   DO i=1,RESNUM2
	   sum1 = 0.0
		DO j=1,RESNUM2
			sum1 = sum1 + delRperbDistMat_XZ(j,i)
		ENDDO
		write(44,*) i, sum1	
		allSumColumn_XZ(i)=sum1/(RESNUM2-1)
	   ENDDO


C############################################################

	   DO J=1,RESNUM
	     DO I=1,52
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=54,61
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=63,66
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=68,96
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=99,108
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=110,146
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=148,149
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=151,186
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=188,208
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=210,261
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=263,270
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=272,275
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=277,305
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=308,317
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=319,355
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=357,358
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=360,395
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=397,417
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=419,470
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=472,479
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=481,484
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=486,514
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=517,526
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=528,564
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=566,567
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=569,604
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=606,626
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=628,679
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=681,688
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=690,693
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=695,723
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=726,735
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=737,773
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=775,776
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=778,813
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	     DO I=815,835
			delRperbDistMat_XZ(i,j) = 0.0
	  	 ENDDO
	   ENDDO	
	   


C############################################################

	   DO i=1,RESNUM2
	   sum = 0.0
		DO j=1,RESNUM2
			sum = sum + delRperbDistMat_XZ(j,i)
		ENDDO
		write(43,*) i, sum	
		bindSumRow_XZ(i) = sum/noPeptideSize*4
	   ENDDO

	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_XZ(j,i)
		ENDDO
	   ENDDO
	
c********************************************************************
c********************************************************************

c                            YZ-direction
	
c********************************************************************
c********************************************************************

	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

    
	DO 1005 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-1,3*nperturbRes
C			delForce(i) = DELF*(ran1(idum)*2-1)
C			delForceRand(i) = delForce(i)
			delForce(i) = 1.0
c		write(*,*) nperturbRes,i,delForce(i)
		ENDDO
		

C******************************************************************
C	Multiply force vector with inverse of BBT = DeltaR
C******************************************************************
c	 DO J=1,RES3
c		write(99,*) j,delForce(j)
c	 ENDDO
	 
	DO I=1,RES3
	 delRperb_YZ(I)= 0.0
	 DO J=1,RES3
	    delRperb_YZ(I)= delRperb_YZ(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
	    XI_YZ(I) = X(I)
        YI_YZ(I) = Y(I)
        ZI_YZ(I) = Z(I)	
	ENDDO


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

C******************************************************************
C	Superimpose the new structure to the initial,
C******************************************************************

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



C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
	IF(PR.GT.PR0)THEN
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
	ENDIF

       sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_YZ(i) = sqrt( delRperb_YZ(3*I-2) *  delRperb_YZ(3*I-2) +  delRperb_YZ(3*I-1) * delRperb_YZ(3*I-1) +  delRperb_YZ(3*I) *  delRperb_YZ(3*I) )
c		   write(*,*) nperturbRes,i,delRperbDist(i)
		   sum = sum + delRperbDist_YZ(i) 		   
	   ENDDO
c	   write(*,*) sum

	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_YZ(i,j) = delRperbDist_YZ(j) / sum
		ENDDO
	   ENDDO

1005	enddo
	
	   DO I=1,RESNUM
		delRperbDistMat_YZ(i,i) = 0.0
	   ENDDO


	   DO I=1,RESNUM
		delRperbDistMat_YZ(i,i+1) = 0.0
		delRperbDistMat_YZ(i,i+2) = 0.0
		delRperbDistMat_YZ(i,i+3) = 0.0
		delRperbDistMat_YZ(i+1,i) = 0.0
		delRperbDistMat_YZ(i+2,i) = 0.0
		delRperbDistMat_YZ(i+3,i) = 0.0
	   ENDDO



	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_YZ(i,j)
		ENDDO
	   ENDDO

	   DO i=1,RESNUM2
	   sum1 = 0.0
		DO j=1,RESNUM2
			sum1 = sum1 + delRperbDistMat_YZ(j,i)
		ENDDO
		write(48,*) i, sum1	
		allSumColumn_YZ(i)=sum1/(RESNUM2-1)
	   ENDDO



C############################################################

	   DO J=1,RESNUM
	     DO I=1,52
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=54,61
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=63,66
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=68,96
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=99,108
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=110,146
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=148,149
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=151,186
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=188,208
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=210,261
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=263,270
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=272,275
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=277,305
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=308,317
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=319,355
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=357,358
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=360,395
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=397,417
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=419,470
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=472,479
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=481,484
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=486,514
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=517,526
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=528,564
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=566,567
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=569,604
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=606,626
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=628,679
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=681,688
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=690,693
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=695,723
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=726,735
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=737,773
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=775,776
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=778,813
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	     DO I=815,835
			delRperbDistMat_YZ(i,j) = 0.0
	  	 ENDDO
	   ENDDO	
	   

C############################################################

	   DO i=1,RESNUM2
	   sum = 0.0
		DO j=1,RESNUM2
			sum = sum + delRperbDistMat_YZ(j,i)
		ENDDO
		write(47,*) i, sum	
		bindSumRow_YZ(i) = sum/noPeptideSize*4
	   ENDDO

	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_YZ(j,i)
		ENDDO
	   ENDDO
	
c********************************************************************
c********************************************************************

c                            XYZ-direction
	
c********************************************************************
c********************************************************************


	DO i=1,NR
		PRCORR(i) = 0.
	ENDDO

   
	DO 1006 nperturbRes = iperturbRes,fperturbRes

		DO i=1,RES3
			delForce(i) = 0.
		ENDDO
	
C	INTERCHANGE THE FOLLOWING IF A NON-RANDOM PERTURBATION DIRECTION IS DESIRED
        DO i=3*nperturbRes-2,3*nperturbRes
C			delForce(i) = DELF*(ran1(idum)*2-1)
C			delForceRand(i) = delForce(i)
			delForce(i) = 1.0
c		write(*,*) nperturbRes,i,delForce(i)
		ENDDO
		

C******************************************************************
C	Multiply force vector with inverse of BBT = DeltaR
C******************************************************************
c	 DO J=1,RES3
c		write(99,*) j,delForce(j)
c	 ENDDO
	 
	DO I=1,RES3
	 delRperb_XYZ(I)= 0.0
	 DO J=1,RES3
	    delRperb_XYZ(I)= delRperb_XYZ(I) + INVBKBT(I,J) * delForce(J)
	  ENDDO
	ENDDO
     
	DO I=1,RESNUM
		XI_XYZ(I) = X(I)
        YI_XYZ(I) = Y(I)
        ZI_XYZ(I) = Z(I)	
	ENDDO


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

C******************************************************************
C	Superimpose the new structure to the initial,
C******************************************************************

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


C	STORE THE STRUCTURE WITH THE BEST PEARSON CORRELATION SCORE
c	INTERCHANGE BETWEEN THE FOLLOWING FOR OTHER CRITERIA OF "BEST"
	IF(PR.GT.PR0)THEN
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
	ENDIF

       sum = 0.0
	   DO i=1,RESNUM
		   delRperbDist_XYZ(i) = sqrt( delRperb_XYZ(3*I-2) *  delRperb_XYZ(3*I-2) +  delRperb_XYZ(3*I-1) * delRperb_XYZ(3*I-1) +  delRperb_XYZ(3*I) *  delRperb_XYZ(3*I) )
c		   write(*,*) nperturbRes,i,delRperbDist(i)
		   sum = sum + delRperbDist_XYZ(i) 		   
	   ENDDO
c	   write(*,*) sum

	   
	   DO I=nperturbRes,nperturbRes
	    DO J=1,RESNUM
			delRperbDistMat_XYZ(i,j) = delRperbDist_XYZ(j) / sum
		ENDDO
	   ENDDO


1006	enddo
	
	   DO I=1,RESNUM
		delRperbDistMat_XYZ(i,i) = 0.0
	   ENDDO


	   DO I=1,RESNUM
		delRperbDistMat_XYZ(i,i+1) = 0.0
		delRperbDistMat_XYZ(i,i+2) = 0.0
		delRperbDistMat_XYZ(i,i+3) = 0.0
		delRperbDistMat_XYZ(i+1,i) = 0.0
		delRperbDistMat_XYZ(i+2,i) = 0.0
		delRperbDistMat_XYZ(i+3,i) = 0.0
	   ENDDO



	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_XYZ(i,j)
		ENDDO
	   ENDDO

	   DO i=1,RESNUM2
	   sum1 = 0.0
		DO j=1,RESNUM2
			sum1 = sum1 + delRperbDistMat_XYZ(j,i)
		ENDDO
		write(55,*) i, sum1	
		allSumColumn_XYZ(i)=sum1/(RESNUM2-1)
	   ENDDO


C############################################################

	   DO J=1,RESNUM
	     DO I=1,52
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=54,61
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=63,66
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=68,96
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=99,108
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=110,146
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=148,149
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=151,186
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=188,208
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=210,261
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=263,270
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=272,275
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=277,305
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=308,317
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=319,355
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=357,358
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=360,395
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=397,417
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=419,470
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=472,479
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=481,484
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=486,514
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=517,526
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=528,564
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=566,567
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=569,604
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=606,626
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=628,679
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=681,688
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=690,693
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=695,723
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=726,735
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=737,773
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=775,776
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=778,813
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	     DO I=815,835
			delRperbDistMat_XYZ(i,j) = 0.0
	  	 ENDDO
	   ENDDO	


C############################################################

	   DO i=1,RESNUM2
	   sum = 0.0
		DO j=1,RESNUM2
			sum = sum + delRperbDistMat_XYZ(j,i)
		ENDDO
		write(54,*) i, sum	
		bindSumRow_XYZ(i) = sum/noPeptideSize*4
	   ENDDO

	   DO i=1,RESNUM
	   sum1 = 0.0
		DO j=1,RESNUM
			sum1 = sum1 + delRperbDistMat_XYZ(j,i)
		ENDDO
	   ENDDO
	

c after oobtained information based on X,Y,Z,XY,XZ,YZ,XYZ directions, we calculate allosteric response ratio for 
c each direction and take the maximum of these values here:
	  write(11,*) 'X,Y,Z,XY,XZ,YZ,XYZ'
	  DO i=1,RESNUM
		ratio_X = (bindSumRow_X(i)/allSumColumn_X(i))
		ratio_Y = (bindSumRow_y(i)/allSumColumn_Y(i))
		ratio_Z = (bindSumRow_z(i)/allSumColumn_Z(i))
		ratio_XY = (bindSumRow_XY(i)/allSumColumn_XY(i))
		ratio_XZ = (bindSumRow_XZ(i)/allSumColumn_XZ(i))
		ratio_YZ = (bindSumRow_YZ(i)/allSumColumn_YZ(i))
		ratio_XYZ = (bindSumRow_XYZ(i)/allSumColumn_XYZ(i))
		write(11,122) i,ratio_X,ratio_Y,ratio_Z,ratio_XY,ratio_XZ,ratio_YZ,ratio_XYZ
		write(12,122) i,MAX(ratio_X,ratio_Y,ratio_Z,ratio_XY,ratio_XZ,ratio_YZ,ratio_XYZ)
	  ENDDO
122   format(i4,7f11.6)



	WRITE(*,*) 'Program finished successfully!'

c666	pause

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

