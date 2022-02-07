
      SUBROUTINE DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                  YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
C     ---------------------------------------------------------
C
      DATA BIG,SMALL /1.0E+20,1.0E-20/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Obtain some viewport and window information about the current 
C    PGPLOT device, without directly accessing the common blocks in
C    pgplot.inc.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    XPERIN   R*4    O       -      Plot X scale in dots/inch.
C    YPERIN   R*4    O       -      Plot Y scale in dots/inch.
C    XOFF     R*4    O       -      Absolute coord of blc of viewport.
C    YOFF     R*4    O       -      Absolute coord of blc of viewport.
C    XLEN     R*4    O       -      Width of viewport in absolute coord.
C    YLEN     R*4    O       -      Height of viewport in absolute coord.
C    XORG     R*4    O       -      Absolute coord of world X=0.
C    YORG     R*4    O       -      Absolute coord of world Y=0.
C    XSCALE   R*4    O       -      Absolute units per world coord in X.
C    YSCALE   R*4    O       -      Absolute units per world coord in Y.
C    XBLC     R*4    O       -      World X coord at blc of window.
C    XTRC     R*4    O       -      World X coord at trc of window.
C    YBLC     R*4    O       -      World Y coord at blc of window.
C    YTRC     R*4    O       -      World Y coord at trc of window.
C
C Globals
C     None.
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGQVP      Inquires about viewport dimensions.
C     PGQWIN     Inquires about world coords of window.
C
C History
C   D. S. Sivia       1 Aug 1996  Initial release.
C-----------------------------------------------------------------------
C
      CALL PGQWIN(XBLC,XTRC,YBLC,YTRC)
      CALL PGQVP(1,XI1,XI2,YI1,YI2)
      CALL PGQVP(3,XOFF,XP2,YOFF,YP2)
      XLEN=ABS(XP2-XOFF)
      YLEN=ABS(YP2-YOFF)
      XPERIN=XLEN/(ABS(XI2-XI1)+SMALL)
      YPERIN=YLEN/(ABS(YI2-YI1)+SMALL)
      XWDIF=XTRC-XBLC
      YWDIF=YTRC-YBLC
      AXWDIF=BIG
      AYWDIF=BIG
      IF (ABS(XWDIF).GT.SMALL) AXWDIF=1.0/XWDIF
      IF (ABS(YWDIF).GT.SMALL) AYWDIF=1.0/YWDIF
      XSCALE=XLEN*AXWDIF
      YSCALE=YLEN*AYWDIF
      XORG=(XOFF*XTRC-XP2*XBLC)*AXWDIF
      YORG=(YOFF*YTRC-YP2*YBLC)*AYWDIF
      END SUBROUTINE DSQINF

      SUBROUTINE PGCELL(A,IDIM,JDIM,I1,I2,J1,J2,FG,BG,TR,NCOLS,R,G,B)
C     ---------------------------------------------------------------
C
      use grpckg1inc
      REAL         A(IDIM,JDIM),TR(6)
      REAL         R(0:NCOLS-1),G(0:NCOLS-1),B(0:NCOLS-1)
      INTEGER      IDIM,JDIM,I1,I2,J1,J2,NCOLS
C
      INTEGER      IR(0:255),IG(0:255),IB(0:255)
      LOGICAL      LPS,LCOLOR,PGNOTO
      CHARACTER*16 TYPE,CHR
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine is designed to do the job of CELL_ARRAY in GKS;
C   that is, it shades elements of a rectangular array with the 
C   appropriate colours passed down in the RGB colour-table. Essentially,
C   it is a colour version of PGGRAY. 
C      The colour-index used for particular array pixel is given by:
C         Colour Index = NINT{[A(i,j)-BG/(FG-BG)]*FLOAT(NCOLS-1)} ,
C   with truncation at 0 and NCOLS-1, as necessary.
C      The transform matrix TR is used to calculate the (bottom left) 
C   world coordinates of the cell which represents each array element:
C         X = TR(1) + TR(2)*I + TR(3)*J
C         Y = TR(4) + TR(5)*I + TR(6)*J  .
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    A        R*4    I   IDIMxJDIM  The array to be plotted.
C    IDIM     I*4    I       -      The first dimension of array A.
C    JDIM     I*4    I       -      The second dimension of array A.
C    I1,I2    I*4    I       -      The inclusive range of the first
C                                   index (I) to be plotted.
C    J1,J2    I*4    I       -      The inclusive range of the second
C                                   index (J) to be plotted.
C    FG       R*4    I       -      The array value which is to appear
C                                   with shade 1 ("foreground").
C    BG       R*4    I       -      The array value which is to appear
C                                   with shade 0 ("background").
C    TR       R*4    I       6      Transformation matrix between array
C                                   grid and world coordinates.
C    NCOLS    I*4    I       -      Number of colours in colour-table.
C    R        R*4    I      NCOLS   Red intensity for colour-table.
C    G        R*4    I      NCOLS   Green intensity for colour-table.
C    B        R*4    I      NCOLS   Blue intensity for colour-table.
C
C Globals
C    GRPCKG1.INC
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGNOTO     Logical function to test if a PGPLOT device is open.
C     GRBPIC     Sends a "begin picture" command to the device driver.
C     PGQCOL     Inquires about the colour capability.
C     PGQCI      Inquires about the current colour index.
C     PGQINF     Inquires about general PGPLOT information.
C     PGBBUF     Recommended initial call (to start a PGPLOT buffer).
C     PGEBUF     Recommended final call (to end a PGPLOT buffer).
C     PGCLPX     Pixel-device support subroutine for PGCELL.
C     PGCLPS     PostScript support subroutine for PGCELL.
C     GREXEC     Dispatches command to appropriate device driver.
C     DSQINF     Inquires about viewport and window dimensions.
C
C History
C   D. S. Sivia       3 Jul 1992  Initial release.
C   D. S. Sivia       6 Feb 1995  Now uses GRGRAY approach instead of
C                                 PGPOLY, and linearly interpolates.
C   D. S. Sivia       6 Mar 1995  Slight changes for Postscript output.
C   D. S. Sivia       1 Aug 1996  Replaced pgplot.inc with DSQINF!
C   D. S. Sivia      21 Oct 1997  Made slightly friendlier for NT.
C   D. S. Sivia      16 Jul 1999  Added a couple of PGPLOT calls to 
C                                 force proper initialisation.
C-----------------------------------------------------------------------
C
C A PGPLOT initialisation precaution.
C
      IF (PGNOTO('PGCELL')) RETURN
      IF (.NOT.GRPLTD(GRCIDE)) CALL GRBPIC
C
C Find out device-type. If not Postscript, then (i) return if less than
C 16 shades available; (ii) save initial colour index (and hope it's 
C less than ICLOW).
C
      NC=256
      LPS=.TRUE.
      CALL PGQINF('TYPE',TYPE,LCHR)
      IF (TYPE.EQ.'PS' .OR. TYPE.EQ.'VPS') THEN
        LCOLOR=.FALSE.
      ELSEIF (TYPE.EQ.'CPS' .OR. TYPE.EQ.'VCPS') THEN
        LCOLOR=.TRUE.
      ELSE
        LPS=.FALSE.
        CALL PGQCOL(IC1,IC2)
        NC=IC2-IC1+1
        IF (NC.LT.16) THEN
          WRITE(*,*)' *** Not enough colours available on this device!'
          RETURN
        ELSE
          ICLOW=4
          IF (NC.GE.96) ICLOW=16
          IC1=IC1+ICLOW
          NC=NC-ICLOW
        ENDIF
        CALL PGQCI(ICSAVE)
        IF (ICSAVE.GE.IC1) ICSAVE=IC1+1
      ENDIF
      CALL PGBBUF
C
C Activate the colour table. If NCOLS is less than the number of colours
C available, simply assign NCOLS; otherwise, use a linear interpolation.
C
      IF (NCOLS.LE.NC) THEN
        NC=NCOLS
        DO 10 I=0,NCOLS-1
          CALL PGSCR(IC1+I,R(I),G(I),B(I))
          IR(I)=NINT(R(I)*255.0)
          IG(I)=NINT(G(I)*255.0)
          IB(I)=NINT(B(I)*255.0)
  10    CONTINUE
      ELSE
        COL=0.0
        DCOL=0.999*FLOAT(NCOLS-1)/FLOAT(NC-1)
        DO 20 I=0,NC-1
          ICOL=INT(COL)
          DICOL=COL-FLOAT(ICOL)
          RL=R(ICOL)+DICOL*(R(ICOL+1)-R(ICOL))
          GL=G(ICOL)+DICOL*(G(ICOL+1)-G(ICOL))
          BL=B(ICOL)+DICOL*(B(ICOL+1)-B(ICOL))
          CALL PGSCR(IC1+I,RL,GL,BL)
          IR(I)=NINT(RL*255.0)
          IG(I)=NINT(GL*255.0)
          IB(I)=NINT(BL*255.0)
          COL=COL+DCOL
  20    CONTINUE
      ENDIF
      ASCALE=FLOAT(NC-1)/(FG-BG)
C
C Check to see whether a pixel device or Postscript is being used, and 
C call the appropriate PGCELL support subrotuine.
C
      IF (LPS) THEN
        CALL PGCLPS(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IR,IG,IB,NC,
     *              LCOLOR)
      ELSE
        NBUF=0
        LCHR=LEN(CHR)
        CALL GREXEC(GRGTYP,4,RBUF,NBUF,CHR,LCHR)
        IF (CHR(7:7).EQ.'P') THEN
          CALL PGCLPX(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IC1,NC)
        ELSE
          WRITE(*,*)' Sorry, PGCELL does not support this device!'
        ENDIF
      ENDIF
C
C Reset the initial colour index.
C
      IF (.NOT. LPS) CALL PGSCI(ICSAVE)
      CALL PGEBUF
      END SUBROUTINE PGCELL

      SUBROUTINE PGCLPX(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IC1,NC)
C     --------------------------------------------------------------
C
C Light-up the device pixels, with colours determined by a linearly
C interpolating array A.
C
      use grpckg1inc
      REAL      A(IDIM,*),TR(*),BUFFER(2050)
      CHARACTER CHR*16
C
      CALL DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *            YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
      IX1=NINT(XOFF)
      IX2=NINT(XOFF+XLEN)
      JY1=NINT(YOFF)
      JY2=NINT(YOFF+YLEN)
      DET=TR(2)*TR(6)-TR(3)*TR(5)
      TR11=+TR(6)/DET
      TR12=-TR(3)/DET
      TR21=-TR(5)/DET
      TR22=+TR(2)/DET
      DX=1.0/XSCALE
      DY=1.0/YSCALE
      X1=(XOFF-XORG)*DX-TR(1)
      Y1=(YOFF-YORG)*DY-TR(4)
      DXX=TR11*DX
      DXY=TR12*DY
      DYX=TR21*DX
      DYY=TR22*DY
      XI1=TR11*X1+TR12*Y1
      YJ1=TR21*X1+TR22*Y1
      DO 40 JY=JY1,JY2
        XI=XI1
        YJ=YJ1
        BUFFER(2)=FLOAT(JY)
        NPIX=2
        DO 30 IX=IX1,IX2
          I=INT(XI)
          J=INT(YJ)
          IF (I.GE.I1 .AND. I.LT.I2 .AND. J.GE.J1 .AND. J.LT.J2) THEN
            IF (NPIX.EQ.2) BUFFER(1)=FLOAT(IX)
            II=I+1
            JJ=J+1
            X=XI-FLOAT(I)
            Y=YJ-FLOAT(J)
            AXY=(1.0-X)*(A(I,J)+Y*(A(I,JJ)-A(I,J)))+
     *               X*(A(II,J)+Y*(A(II,JJ)-A(II,J)))
            K=NINT((AXY-BG)*ASCALE)
            IF (K.LT.0) K=0
            IF (K.GE.NC) K=NC-1
            NPIX=NPIX+1
            BUFFER(NPIX)=FLOAT(K+IC1)
          ENDIF
          XI=XI+DXX
          YJ=YJ+DYX
  30    CONTINUE
        CALL GREXEC(GRGTYP,26,BUFFER,NPIX,CHR,LCHR)
        XI1=XI1+DXY
        YJ1=YJ1+DYY
  40  CONTINUE
      END SUBROUTINE PGCLPX

      SUBROUTINE PGCLPS(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IR,IG,IB,
     *                  NC,LCOLOR)
C     -----------------------------------------------------------------
C
C Postscript support subroutine for PGCELL.
C
      use grpckg1inc
      REAL      A(IDIM,*),TR(*)
      INTEGER   IR(0:*),IG(0:*),IB(0:*)
      INTEGER   VALUE(33)
      LOGICAL   LCOLOR
      CHARACTER INLINE*80
C
C Set clipping rectangle in device.
C
      CALL DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *            YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
      WRITE (INLINE,100) NINT(XOFF),NINT(YOFF),NINT(XLEN),NINT(YLEN),
     *          -NINT(XLEN)
 100  FORMAT(I6,I6,' moveto ',I6, ' 0 rlineto  0 ',I6,' rlineto ',
     *I6,' 0 rlineto')
      CALL GRTERM
      CALL GRESC(' newpath ')
      CALL GRESC(INLINE)
      CALL GRESC(' closepath ')
C
C Work out the nunmber of X and Y pixels for PS image, with NDOTS per 
C inch, and build an image transformation matrix.
C
      NDOTS=100
C      IF (LCOLOR) NDOTS=50
      NXP=NINT(FLOAT(NDOTS)*XLEN/XPERIN)
      NYP=NINT(FLOAT(NDOTS)*YLEN/YPERIN)
      DXI=FLOAT(I2-I1)/FLOAT(NXP)
      DYJ=FLOAT(J2-J1)/FLOAT(NYP)
      DET=TR(2)*TR(6)-TR(3)*TR(5)
      TR11=+TR(6)/DET
      TR12=-TR(3)/DET
      TR21=-TR(5)/DET
      TR22=+TR(2)/DET
      AT=TR11/(XSCALE*DXI)
      BT=TR21/(XSCALE*DYJ)
      CT=TR12/(YSCALE*DXI)
      DT=TR22/(YSCALE*DYJ)
      TX=-(TR11*(XORG/XSCALE+TR(1))+TR12*(YORG/YSCALE+TR(4))+I1)/DXI
      TY=-(TR21*(XORG/XSCALE+TR(1))+TR22*(YORG/YSCALE+TR(4))+J1)/DYJ
C
C Use a PostScript "image" operator.
C
      WRITE (INLINE, '(A,I5,A)') '/picstr ',NXP,' string def'
      CALL GRESC(INLINE)
      WRITE (INLINE,110) NXP,NYP,AT,BT,CT,DT,TX,TY
 110  FORMAT(2I4,' 8 [',6(1PE10.3,' '),']')
      CALL GRESC(INLINE)
      CALL GRESC('{ currentfile picstr readhexstring pop}')
      IF (LCOLOR) THEN
        CALL GRESC('  false 3 colorimage')
      ELSE
        CALL GRESC('  image')
      ENDIF
C
C Write out the image array in hexadecimal.
C
      ASCALE=ASCALE*255.0/FLOAT(NC-1)
      YJ=FLOAT(J1)
      DO 20 JP=1,NYP
        J=INT(YJ)
        Y=YJ-FLOAT(J)
        JJ=J+1
        IF (JJ.GT.J2) JJ=J2
        XI=FLOAT(I1)
        L=0
        DO 10 IP=1,NXP
          L=L+1
          I=INT(XI)
          X=XI-FLOAT(I)
          II=I+1
          IF (II.GT.I2) II=I2
          AXY=(1.0-X)*(A(I,J)+Y*(A(I,JJ)-A(I,J)))+
     *             X*(A(II,J)+Y*(A(II,JJ)-A(II,J)))
          IC=NINT((AXY-BG)*ASCALE)
          IF (IC.LT.0) IC=0
          IF (IC.GE.NC) IC=NC-1
          IF (LCOLOR) THEN
            VALUE(L)=IR(IC)
            VALUE(L+1)=IG(IC)
            VALUE(L+2)=IB(IC)
            L=L+2
          ELSE
            VALUE(L)=(IR(IC)+IG(IC)+IB(IC))/3
          ENDIF
          IF (L.EQ.33) THEN
            WRITE(INLINE,120) (VALUE(K),K=1,33)
 120        FORMAT(33Z2.2)
            CALL GRESC(INLINE(1:66))
            L=0
          ENDIF
          XI=XI+DXI
   10   CONTINUE
        IF (L.NE.0) THEN
          WRITE(INLINE,120) (VALUE(K),K=1,L)
          CALL GRESC(INLINE(1:2*L))
        ENDIF
        YJ=YJ+DYJ
   20 CONTINUE
      CALL GRESC(' newpath ')
      CALL GRTERM
      END SUBROUTINE PGCLPS

      SUBROUTINE PLOT(X,NX,Y,NY,Z,N1,N2,W,SIZE,IWIDTH,XLBL,YLBL,TITL)
C     ---------------------------------------------------------------
C
      use coltabs
      REAL          X(*),Y(*),Z(N1,N2)
      CHARACTER*(*) XLBL,YLBL,TITL
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C Purpose
C      This subroutine plots "data" defined on a regularly-spaced 
C   rectangular grid of points Z(I,J). With the default choice for the 
C   PGCELL routine that is linked, the output is a linearly-interpolated 
C   map (rather than coarse rectangular boxes).
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    Z        R*4    I    N1 x N2   The recatangular "data"-array.
C    Z        R*4    O    N1 x N2   A scaled, and clipped, version of 
C                                   the input array(!).
C    N1       I*4    I       -      The first dimension of array Z.
C    N2       I*4    I       -      The second dimension of array Z.
C    X        R*4    I      NX      Array of X-coordinates.
C    NX       I*4    I       -      Number of X-pixels to be plotted 
C                                  (usually = N1, but must be <= N1).
C    Y        R*4    I      NY      Array of Y-coordinates.
C    NY       I*4    I       -      Number of Y-pixels to be plotted 
C                                  (usually = N2, but must be <= N2).
C    W        R*4    I       -      Dummy parameter (back-compatibility).
C    SIZE     R*4    I       -      Character-size for plot (try 1.5).
C    IWIDTH   I*4    I       -      Line-width for plot (try 2).
C    XLBL     A*1    I     *(*)     Label for X-axis.
C    YLBL     A*1    I     *(*)     Label for Y-axis.
C    TITL     A*1    I     *(*)     Title for plot.
C
C Globals
C    COLTABS.INC
C
C History
C   Initial release.                                    DSS:  3 Jul 1992
C   Minor changes to conform with new PGCELL.           DSS:  6 Feb 1995
C   Put in option to over-lay contours.                 DSS: 21 Feb 1995
C   Now has proper 3-d surface surface rendering.       DSS: 27 Aug 1997
C   Fortran made LINUX-friendly!                        DSS: 15 Sep 1997
C   Choose white background colour for postscript.      DSS:  5 Feb 1999
C-----------------------------------------------------------------------
C
C      INCLUDE  'COLTABS.INC'
      REAL      RGB(3,3),EYE(3),LIGHT(3),LATICE(3,3),LUTUSR(3,256)
      CHARACTER STRING*32,TYPE*16,CHR*16
      LOGICAL   OVRLAY,LSHIN
      DATA      EYE,LIGHT /0.0,0.0,1000.0,-1.0,-1.0,-1.0/
      DATA      RGB /0.0,0.0,1.0,0.35,0.35,0.35,1.0,1.0,1.0/
C
      WRITE(*,*)
      WRITE(*,*)'                (0) Contour'
      WRITE(*,*)'                (1) Surface'
      WRITE(*,*)'                (2) Colour: Grey-Scale'
      WRITE(*,*)'                (3) Colour: Heat'
      WRITE(*,*)'                (4) Colour: Rainbow Spectrum'
      WRITE(*,*)'                (5) Colour: BGYRW'
      WRITE(*,*)'                (6) Colour: Serpent'
      WRITE(*,*)'                (7) Colour: Read in from file'
      WRITE(*,*)
   1  IPLOT=0
      WRITE(*,100)
 100  FORMAT(' PLOT>  Type ?  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) IPLOT
      IF (IPLOT.EQ.0) THEN
        CALL CONTOR(Z,NX,NY,X,Y,N1,N2,SIZE,IWIDTH,.FALSE.)
      ELSEIF (IPLOT.EQ.1) THEN
        IF (N1.NE.NX .OR. N2.NE.NY) THEN
          WRITE(*,*)' Sorry folks, SURFACE option needs N1=NX & N2=NY!'
          WRITE(*,*)
          GOTO 1
        ENDIF
        CALL SRFCOL(RGB,NCB,ICTAB,DIFUS,SHIN,POLISH,LSHIN,3,7)
        CALL EULER(LATICE)
        CALL FMXMN(Z,N1*N2,DHIGH,DLOW,DOFSET)
        CALL PGBEGIN(0,'?',1,1)
        CALL PGPAPER(0.0,1.0)
        CALL PGQCOL(ICMIN,ICMAX)
        NCBAND=MIN(NCB,ICMAX-17+1)
        CALL PGSCH(SIZE)
        CALL PGSLW(IWIDTH)
        CALL PGVPORT(0.0,1.0,0.0,1.0)
        CALL PGWINDOW(-0.87,0.92,-0.87,0.87)
        CALL PGSCI(0)
        CALL PGBOX('BC',0.0,0,'BC',0.0,0)
        CALL PGSCI(1)
        CALL PGQINF('TYPE',TYPE,LCHR)
        IF (TYPE.EQ.'PS' .OR. TYPE.EQ.'VPS' .OR. TYPE.EQ.'CPS' .OR.
     *      TYPE.EQ.'VCPS') THEN
          CALL SBFINT(RGB(1,3),16,1,1,MAXBUF)
        ELSE
          CALL SBFINT(RGB(1,2),16,1,1,MAXBUF)
        ENDIF
        IF (ICTAB.LE.2) THEN
          CALL COLINT(RGB,17,ICMAX,DIFUS,SHIN,POLISH)
        ELSEIF (ICTAB.EQ.3) THEN
          CALL COLSRF(HEAT,256,1.0,17,ICMAX,NCB,DIFUS,SHIN,POLISH)
        ELSEIF (ICTAB.EQ.4) THEN
          CALL COLSRF(SPECTRUM(1,2),255,1.0,17,ICMAX,NCB,DIFUS,SHIN,
     *                POLISH)
        ELSEIF (ICTAB.EQ.5) THEN
          CALL COLSRF(BGYRW,256,1.0,17,ICMAX,NCB,DIFUS,SHIN,POLISH)
        ELSEIF (ICTAB.EQ.6) THEN
          CALL COLSRF(SERP,256,1.0,17,ICMAX,NCB,DIFUS,SHIN,POLISH)
        ELSE
          NCLMAX=256
          CALL LUTIN(LUTUSR,NCLMAX,NCLUSR,IFLAG)
          CALL COLSRF(LUTUSR,NCLUSR,1.0,17,ICMAX,NCB,DIFUS,SHIN,POLISH)
        ENDIF
        CALL SB2SRF(EYE,LATICE,Z,N1-1,N2-1,DLOW,DHIGH,1.0,17,ICMAX,NCB,
     *              LIGHT,LSHIN)
        CALL AXES3D(EYE,LATICE,X(1),X(N1),Y(1),Y(N2),XLBL,YLBL,SIZE,
     *              DLOW,DHIGH,DOFSET,Z(1,1),Z(N1,1),Z(N1,N2),Z(1,N2))
        CALL SBFCLS(1)
        CALL PGMTEXT('T',-1.2,0.5,0.5,TITL)
      ELSEIF (IPLOT.LE.7) THEN
        CALL GREY(Z,NX,NY,X,Y,N1,N2,IPLOT,SIZE,IWIDTH)
      ELSE
        GOTO 1
      ENDIF
      IF (IPLOT.NE.1) CALL PGLABEL(XLBL,YLBL,TITL)
      IF (IPLOT.GT.1) THEN
        WRITE(*,*)
        CALL LOGQYN(' PLOT> Over-lay contours ?','N',OVRLAY)
        IF (OVRLAY) CALL CONTOR(Z,NX,NY,X,Y,N1,N2,SIZE,IWIDTH,OVRLAY)
      ENDIF
      CALL PGEND
      END SUBROUTINE PLOT

      SUBROUTINE FMXMN(Y,N,YMAX,YMIN,YOFSET)
C     --------------------------------------
C
      REAL         Y(*)
      CHARACTER*32 STRING
C
      YMIN1=+1.0E25
      YMAX1=-1.0E25
      DO 10 I=1,N
        YMIN1=MIN(YMIN1,Y(I))
        YMAX1=MAX(YMAX1,Y(I))
  10  CONTINUE
   1  YMIN=YMIN1
      YMAX=YMAX1
      WRITE(*,100) YMIN,YMAX
 100  FORMAT(' Surface>  Zmin & Zmax for plot ? (def=',1pe10.3,
     *       ',',E10.3,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) YMIN,YMAX
      IF (YMAX.LE.YMIN) GOTO 1
      YOFSET=MAX(-YMIN,0.0)
      IF (YOFSET.GT.0.0) THEN
        YMIN=YMIN+YOFSET
        YMAX=YMAX+YOFSET
        DO 20 I=1,N
  20      Y(I)=Y(I)+YOFSET
      ENDIF
      END SUBROUTINE FMXMN

      SUBROUTINE EULER(LATICE)
C     ------------------------
C
      CHARACTER*32 STRING
      REAL         LATICE(3,*)
      DATA         PIRAD /0.01745329252/
C
   1  WRITE(*,100)
 100  FORMAT(' Surface> Rotation and tilt ?  (def=45,30 deg) : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.EQ.0) THEN
        SINA=-0.7071067814
        COSA=+0.7071067814
        SINB=-0.5
        COSB=0.866025404
      ELSE
        READ(STRING,*,ERR=1) IA,IB
        SINA=-SIN(FLOAT(IA)*PIRAD)
        COSA=COS(FLOAT(IA)*PIRAD)
        SINB=-SIN(FLOAT(IB)*PIRAD)
        COSB=COS(FLOAT(IB)*PIRAD)
      ENDIF
      CALL ROTY(-0.5,-0.5,+0.5,U,V,W,SINA,COSA)
      CALL ROTX(U,V,W,LATICE(1,1),LATICE(2,1),Z1,SINB,COSB)
      LATICE(3,1)=Z1-1.0
      CALL ROTY(+0.5,-0.5,+0.5,U,V,W,SINA,COSA)
      CALL ROTX(U,V,W,LATICE(1,2),LATICE(2,2),Z2,SINB,COSB)
      LATICE(3,2)=Z2-1.0
      CALL ROTY(-0.5,-0.5,-0.5,U,V,W,SINA,COSA)
      CALL ROTX(U,V,W,LATICE(1,3),LATICE(2,3),Z3,SINB,COSB)
      LATICE(3,3)=Z3-1.0
      END SUBROUTINE EULER

      SUBROUTINE ROTX(X,Y,Z,U,V,W,S,C)
C     --------------------------------
C
      U=X
      V=Y*C+Z*S
      W=-Y*S+Z*C
      END SUBROUTINE ROTX

      SUBROUTINE ROTY(X,Y,Z,U,V,W,S,C)
C     --------------------------------
C
      U=X*C-Z*S
      V=Y
      W=X*S+Z*C
      END SUBROUTINE ROTY

      SUBROUTINE SRFCOL(RGB,NCBAND,ICTAB,DIF,SHIN,POLISH,LSHIN,IC1,IC2)
C     -----------------------------------------------------------------
C
      REAL         RGB(*)
      LOGICAL      LSHIN
      CHARACTER*32 STRING
C
      POLISH=1.0
   1  ICTAB=0
      WRITE(*,100) IC1,IC2,ICTAB
 100  FORMAT(' Surface> Colour table ?  (',I1,'-',I1,',def=',I1,
     *       ') : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) ICTAB
      IF (ICTAB.GT.IC2) THEN
        GOTO 1
      ELSEIF (ICTAB.LT.IC1) THEN
        NCBAND=1
   2    RGB(1)=0.0
        RGB(2)=1.0
        RGB(3)=0.0
        WRITE(*,110) (RGB(I),I=1,3)
 110    FORMAT(' Surface> RGB colour ?  (def=',F3.1,',',F3.1,',',
     *    F3.1,') : ',$)
        CALL FORMTQ(STRING,32,NNN)
        IF (NNN.NE.0) READ(STRING,*,ERR=2) (RGB(I),I=1,3)
        IF (RGB(1)+RGB(2)+RGB(3).LE.0.05) GOTO 2 
      ELSE
   3    NCBAND=8
        WRITE(*,120) NCBAND
 120    FORMAT(' Surface> No. of colour-bands ?  (def=',I1,') : ',$)
        CALL FORMTQ(STRING,32,NNN)
        IF (NNN.NE.0) READ(STRING,*,ERR=3) NCBAND
        NCBAND=MAX(MIN(NCBAND,64),1)
      ENDIF
      CALL LOGQYN(' Surface> A shiny gloss ?','N',LSHIN)
      IF (LSHIN) THEN
        SHIN=1.0
        DIF=0.0
      ELSE
        SHIN=0.0
   4    DIF=0.7
        WRITE(*,130) DIF
 130    FORMAT(' Surface> Diffusiveness ?  (def=',F3.1,') : ',$)
        CALL FORMTQ(STRING,32,NNN)
        IF (NNN.NE.0) READ(STRING,*,ERR=4) DIF
        DIFUSE=MAX(MIN(DIF,1.0),0.1)
      ENDIF
      END SUBROUTINE SRFCOL

      SUBROUTINE AXES3D(EYE,LATICE,XMIN,XMAX,YMIN,YMAX,XLBL,YLBL,SIZE,
     *                  DLOW,DHIGH,DOFSET,D00,DX0,DXY,D0Y)
C     ----------------------------------------------------------------
C
      REAL          EYE(*),LATICE(3,*)
      REAL          PIVX(3),PIVY(3),ORX(3,2),ORY(3,2),LATCAB(3,4)
      CHARACTER*(*) XLBL,YLBL
C
      IF (XMAX.LE.XMIN .OR. YMAX.LE.YMIN) RETURN
      SCLA=0.15*SIZE
      AX=LATICE(1,2)-LATICE(1,1)
      AY=LATICE(2,2)-LATICE(2,1)
      AZ=LATICE(3,2)-LATICE(3,1)
      BX=LATICE(1,3)-LATICE(1,1)
      BY=LATICE(2,3)-LATICE(2,1)
      BZ=LATICE(3,3)-LATICE(3,1)
      CX=AY*BZ-BY*AZ
      CY=AZ*BX-BZ*AX
      CZ=AX*BY-BX*AY
      XSIGN=+1.0
      XSCL=-SCLA
      IF (CY*BZ.GT.0.0) THEN
        XSIGN=-1.0
        XSCL=1.0+SCLA
      ENDIF
      ORX(1,1)=XSIGN*AX
      ORX(2,1)=XSIGN*AY
      ORX(3,1)=XSIGN*AZ
      ORX(1,2)=XSIGN*BX
      ORX(2,2)=XSIGN*BY
      ORX(3,2)=XSIGN*BZ
      PIVX(1)=0.5*(LATICE(1,1)+LATICE(1,2))+XSCL*BX
      PIVX(2)=0.5*(LATICE(2,1)+LATICE(2,2))+XSCL*BY
      PIVX(3)=0.5*(LATICE(3,1)+LATICE(3,2))+XSCL*BZ
      CALL SBTEXT(EYE,XLBL,1,PIVX,0.5,ORX,SCLA*0.2)
      CALL AXNUMS(EYE,XMIN,XMAX,PIVX,ORX,SCLA,XSIGN)
      YSIGN=-1.0
      YSCL=-SCLA
      IF (CY*AZ.GT.0.0) THEN
        YSIGN=1.0
        YSCL=1.0+SCLA
      ENDIF
      ORY(1,1)=YSIGN*BX
      ORY(2,1)=YSIGN*BY
      ORY(3,1)=YSIGN*BZ
      ORY(1,2)=-YSIGN*AX
      ORY(2,2)=-YSIGN*AY
      ORY(3,2)=-YSIGN*AZ
      PIVY(1)=0.5*(LATICE(1,1)+LATICE(1,3))+YSCL*AX
      PIVY(2)=0.5*(LATICE(2,1)+LATICE(2,3))+YSCL*AY
      PIVY(3)=0.5*(LATICE(3,1)+LATICE(3,3))+YSCL*AZ
      CALL SBTEXT(EYE,YLBL,1,PIVY,0.5,ORY,SCLA*0.2)
      CALL AXNUMS(EYE,YMIN,YMAX,PIVY,ORY,SCLA,YSIGN)
      LATCAB(1,1)=LATICE(1,2)+BX
      LATCAB(2,1)=LATICE(2,2)+BY
      LATCAB(3,1)=LATICE(3,2)+BZ
      CALL SBLINE(EYE,LATICE(1,1),LATICE(1,2),1,.FALSE.)
      CALL SBLINE(EYE,LATICE(1,2),LATCAB(1,1),1,.FALSE.)
      CALL SBLINE(EYE,LATCAB(1,1),LATICE(1,3),1,.FALSE.)
      CALL SBLINE(EYE,LATICE(1,3),LATICE(1,1),1,.FALSE.)
      ZSCALE=1.0/MAX(DHIGH-DLOW,1.0E-20)
      FRACZ=MAX((D00-DLOW)*ZSCALE,0.0)
      LATCAB(1,2)=LATICE(1,1)+FRACZ*CX
      LATCAB(2,2)=LATICE(2,1)+FRACZ*CY
      LATCAB(3,2)=LATICE(3,1)+FRACZ*CZ
      CALL SBLINE(EYE,LATICE(1,1),LATCAB(1,2),1,.FALSE.)
      CALL VCOPY(LATICE(1,1),LATCAB(1,3),3)
      ZXMIN=LATCAB(1,3)
      FRACZ=MAX((DX0-DLOW)*ZSCALE,0.0)
      LATCAB(1,2)=LATICE(1,2)+FRACZ*CX
      LATCAB(2,2)=LATICE(2,2)+FRACZ*CY
      LATCAB(3,2)=LATICE(3,2)+FRACZ*CZ
      CALL SBLINE(EYE,LATICE(1,2),LATCAB(1,2),1,.FALSE.)
      IF (LATICE(1,2).GT.ZXMIN) THEN
        CALL VCOPY(LATICE(1,2),LATCAB(1,3),3)
        ZXMIN=LATCAB(1,3)
      ENDIF
      FRACZ=MAX((DXY-DLOW)*ZSCALE,0.0)
      LATCAB(1,2)=LATCAB(1,1)+FRACZ*CX
      LATCAB(2,2)=LATCAB(2,1)+FRACZ*CY
      LATCAB(3,2)=LATCAB(3,1)+FRACZ*CZ
      IF (LATCAB(1,1).GT.ZXMIN) THEN
        CALL VCOPY(LATCAB(1,1),LATCAB(1,3),3)
        ZXMIN=LATCAB(1,3)
      ENDIF
      CALL SBLINE(EYE,LATCAB(1,1),LATCAB(1,2),1,.FALSE.)
      FRACZ=MAX((D0Y-DLOW)*ZSCALE,0.0)
      LATCAB(1,2)=LATICE(1,3)+FRACZ*CX
      LATCAB(2,2)=LATICE(2,3)+FRACZ*CY
      LATCAB(3,2)=LATICE(3,3)+FRACZ*CZ
      CALL SBLINE(EYE,LATICE(1,3),LATCAB(1,2),1,.FALSE.)
      IF (LATICE(1,3).GT.ZXMIN) THEN
        CALL VCOPY(LATICE(1,3),LATCAB(1,3),3)
        ZXMIN=LATCAB(1,3)
      ENDIF
      LATCAB(1,4)=LATCAB(1,3)+CX
      LATCAB(2,4)=LATCAB(2,3)+CY
      LATCAB(3,4)=LATCAB(3,3)+CZ
      CALL SBLINE(EYE,LATCAB(1,3),LATCAB(1,4),1,.FALSE.)
      CALL AZNUMS(EYE,DLOW-DOFSET,DHIGH-DOFSET,LATCAB(1,3),SCLA)
      END SUBROUTINE AXES3D

      SUBROUTINE AXNUMS(EYE,XMIN,XMAX,PIVX,ORX,SCLA,XSIGN)
C     ----------------------------------------------------
C
      REAL      EYE(*),PIVX(*),ORX(3,*)
      REAL      END1(3),END2(3),PIVOT(3)
      CHARACTER NLBL*20
      DATA      FRTICK,FRNUM /0.02,0.10/
C
      XR=PGRND(XMAX-XMIN,NSUB)
      DX=XR/FLOAT(NSUB)
      IF (DX.LE.1.0E-20) RETURN
   1  XJ=DX*FLOAT(1+INT(XMIN/DX))
      IF ((XJ+DX).GE.XMAX) THEN
        DX=DX/2.0
        NSUB=NSUB*2
        GOTO 1
      ENDIF
      IF (XMIN.LT.0.0) XJ=XJ-DX
      XN=XSIGN/(XMAX-XMIN)
      XH=0.5*(XMIN+XMAX)
      DO 20 J=1,NSUB
        IF (XJ.GT.XMAX) RETURN
        XF=XN*(XJ-XH)
        DO 10 I=1,3
          END1(I)=PIVX(I)+XF*ORX(I,1)+SCLA*ORX(I,2)
          END2(I)=END1(I)-FRTICK*ORX(I,2)
          PIVOT(I)=END1(I)-FRNUM*ORX(I,2)
  10    CONTINUE
        IPOWER=INT(LOG10(ABS(XJ)+1.0E-10))-5
        IF (XJ.LT.1.0) IPOWER=IPOWER-1
        X=XJ/(10.0**IPOWER)
        IMANTS=NINT(X)
        CALL PGNUMB(IMANTS,IPOWER,0,NLBL,NC)
        CALL SBLINE(EYE,END1,END2,1,.FALSE.)
        CALL SBTEXT(EYE,NLBL,1,PIVOT,0.5,ORX,SCLA*0.15)
        XJ=XJ+DX
  20  CONTINUE
      END SUBROUTINE AXNUMS

      SUBROUTINE AZNUMS(EYE,ZMIN,ZMAX,LATZ,SCLA)
C     ------------------------------------------
C
      REAL      EYE(*),LATZ(3,*)
      REAL      END1(3),END2(3),PIVOT(3),ORIENT(3,3)
      CHARACTER NLBL*20
      DATA      FRTICK,FRNUM /0.05,0.10/
C
      ORIENT(1,1)=SCLA
      ORIENT(2,1)=0.0
      ORIENT(3,1)=0.0
      ZSIGN=+1.0
      IF ((LATZ(2,2)-LATZ(2,1)).LT.0.0) ZSIGN=-1.0
      DO 10 I=1,3
        ORIENT(I,3)=LATZ(I,2)-LATZ(I,1)
        ORIENT(I,2)=ZSIGN*SCLA*ORIENT(I,3)
  10  CONTINUE
      ZR=PGRND(ZMAX-ZMIN,NSUB)
      DZ=ZR/FLOAT(NSUB)
      IF (DZ.LE.1.0E-20) RETURN
   1  ZJ=DZ*FLOAT(1+INT(ZMIN/DZ))
      IF ((ZJ+DZ).GE.ZMAX) THEN
        DZ=DZ/2.0
        NSUB=NSUB*2
        GOTO 1
      ENDIF
      IF (ZMIN.LE.0.0) ZJ=ZJ-DZ
      ZN=1.0/(ZMAX-ZMIN)
      DO 30 J=1,NSUB
        IF (ZJ.GT.ZMAX) RETURN
        ZF=ZN*(ZJ-ZMIN)
        DO 20 I=1,3
          END1(I)=LATZ(I,1)+ZF*ORIENT(I,3)
          END2(I)=END1(I)+FRTICK*ORIENT(I,1)
          PIVOT(I)=END1(I)+FRNUM*ORIENT(I,1)-FRTICK*ORIENT(I,2)
  20    CONTINUE
        IPOWER=INT(LOG10(ABS(ZJ)+1.0E-10))-5
        IF (ZJ.LT.1.0) IPOWER=IPOWER-1
        Z=ZJ/(10.0**IPOWER)
        IMANTS=NINT(Z)
        CALL PGNUMB(IMANTS,IPOWER,0,NLBL,NC)
        CALL SBLINE(EYE,END1,END2,1,.FALSE.)
        CALL SBTEXT(EYE,NLBL,1,PIVOT,0.0,ORIENT,SCLA*0.12)
        ZJ=ZJ+DZ
  30  CONTINUE
      END SUBROUTINE AZNUMS

      SUBROUTINE CONTOR(MAP,NX,NY,X,Y,N1,N2,SIZE,IWIDTH,OVRLAY)
C     ---------------------------------------------------------
C
      REAL    MAP(*),X(*),Y(*)
      REAL    CONT(25,2),TR(6)
      INTEGER NCONT(2),LCOLOR(2),LSTYLE(2),LWIDTH(2)
      LOGICAL OVRLAY
      DATA    LCOLOR,LSTYLE,LWIDTH,NWIND,IDP /1,3,1,1,1,2,1,0/
C
      CALL INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,N1*N2)
      IF (OVRLAY) THEN
        CALL AUTCNT(NCONT,CONT,DMIN,DMAX,OVRLAY)
        CALL PGSLW(1)
        CALL PGCONT(MAP,N1,N2,1,NX,1,NY,CONT,NCONT,TR)
        CALL PGSLW(IWIDTH)
        RETURN
      ENDIF
      CALL ASKCNT(NCONT,CONT,DMIN,DMAX)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGENV(X(1),X(NX),Y(1),Y(NY),0,0)
      DO 20 I=1,2
        CALL PGSCI(LCOLOR(I))
        CALL PGSLS(LSTYLE(I))
        LWDTH=2*LWIDTH(I)
        IF (LWDTH.GT.7) LWDTH=7
        CALL PGSLW(LWDTH)
        CALL PGCONT(MAP,N1,N2,1,NX,1,NY,CONT(1,I),NCONT(I),TR)
  20  CONTINUE
      CALL PGSCI(1)
      CALL PGSLW(IWIDTH)
      END SUBROUTINE CONTOR

      SUBROUTINE INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,NMAP)
C     ------------------------------------------------
C
      REAL X(*),Y(*),TR(*),MAP(*)
C
      TR(1)=X(1)-(X(NX)-X(1))/FLOAT(NX-1)
      TR(2)=(X(NX)-X(1))/FLOAT(NX-1)
      TR(3)=0.0
      TR(4)=Y(1)-(Y(NY)-Y(1))/FLOAT(NY-1)
      TR(5)=0.0
      TR(6)=(Y(NY)-Y(1))/FLOAT(NY-1)
      DMIN=+1.0E+20
      DMAX=-1.0E+20
      DTOT=0.0
      DO 10 I=1,NMAP
        F=MAP(I)
        DTOT=DTOT+F
        IF (F.GT.DMAX) DMAX=F
        IF (F.LT.DMIN) DMIN=F
  10  CONTINUE
      WRITE(*,*)
      WRITE(*,*)' Minimum value = ',DMIN
      WRITE(*,*)' Maximum value = ',DMAX
      WRITE(*,*)' Total flux    = ',DTOT
      WRITE(*,*)
      END SUBROUTINE INIT

      SUBROUTINE ASKCNT(NCONT,CONT,DMIN,DMAX)
C     ---------------------------------------
C
      REAL           CONT(25,*),C(25)
      INTEGER        NCONT (*)
      LOGICAL        AUTO
      CHARACTER*132 STRING
C
      WRITE(*,*)
      CALL LOGQYN(' Contours>  Autoscale (linear) ?','Y',AUTO)
      IF (AUTO) THEN
        CALL AUTCNT(NCONT,CONT,DMIN,DMAX,.FALSE.)
        RETURN
      ENDIF
      WRITE(*,*)
      WRITE(*,*)'          ***** Contour Values  *****'
      WRITE(*,*)
   1  WRITE(*,100)
 100  FORMAT(' Thin>  ',$)
      READ(*,200,ERR=1) STRING
 200  FORMAT(A)
      CALL FINDNC(STRING,132,NCONT(1))
      READ(STRING,*,ERR=1) (CONT(I,1),I=1,NCONT(1))
      WRITE(*,*)
   2  WRITE(*,110)
 110  FORMAT(' Thick>  ',$)
      READ(*,200,ERR=2) STRING
      CALL FINDNC(STRING,132,NCONT(2))
      READ(STRING,*,ERR=2) (CONT(I,2),I=1,NCONT(2))
      END SUBROUTINE ASKCNT

      SUBROUTINE FINDNC(STRING,NCHARS,NC)
C     -----------------------------------
C
      CHARACTER STRING(*)
C
      DO 10 I=1,NCHARS
  10    IF (STRING(I).NE.' ') GOTO 1
   1  IMIN=I
      DO 20 I=NCHARS,1,-1
  20    IF (STRING(I).NE.' ') GOTO 2
   2  IMAX=I
      IF (IMIN.LE.IMAX) THEN
        NC=1
        J=IMIN
        DO 30 I=IMIN,IMAX
          IF (J.GT.IMAX) RETURN
          IF (STRING(J).EQ.' ' .OR. STRING(J).EQ.',') THEN
            NC=NC+1
   3        IF (STRING(J+1).EQ.' ' .OR. STRING(J+1).EQ.',') THEN
              J=J+1
              GOTO 3
            ENDIF
          ENDIF
          J=J+1
  30    CONTINUE
      ELSE
        NC=0
      ENDIF
      END SUBROUTINE FINDNC

      SUBROUTINE AUTCNT(NCONT,CONT,DMIN,DMAX,OVRLAY)
C     ----------------------------------------------
C
      REAL         CONT(25,*)
      INTEGER      NCONT(*)
      LOGICAL       OVRLAY
      CHARACTER*32 STRING
C
      WRITE(*,*)
   1  WRITE(*,100) DMIN,DMAX
 100  FORMAT(' Autoscale> Range ? (def=',1pe10.3,' to ',e10.3,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) DMIN1,DMAX1
      IF (NNN.NE.0) THEN
        DMIN=DMIN1
        DMAX=DMAX1
      ENDIF
   2  N=10
      WRITE(*,110) N
 110  FORMAT(' Autoscale>  No. of contours ?  (def=',I2,') : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=2) N
      N=MAX(MIN(N,50),2)
      NCONT(1)=(N+1)/2
      NCONT(2)=N-NCONT(1)
      IF (OVRLAY) THEN
        IF (N.GT.25) N=25
        NCONT(1)=N
        NCONT(2)=0
      ENDIF
      XINC=(DMAX-DMIN)/FLOAT(N)
      X=DMIN+0.5*XINC
      DO 10 I=1,NCONT(1)
        CONT(I,1)=X
        X=X+XINC
  10  CONTINUE
      DO 20 I=1,NCONT(2)
        CONT(I,2)=X
        X=X+XINC
  20  CONTINUE
      END SUBROUTINE AUTCNT

      SUBROUTINE GREY(MAP,NX,NY,X,Y,N1,N2,IPLOT,SIZE,IWIDTH)
C     ------------------------------------------------------
C
      REAL MAP(*),X(*),Y(*)
      REAL TR(6),TRCOL(6),COLBAR(2,256),LUTUSR(3,256)
      REAL RED(256),GR(256),BL(256)
      DATA TRCOL /0.0,1.0,0.0,0.0,0.0,1.0/
C
      NCOLS=256
      CALL INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,N1*N2)
      CALL STGREY(MAP,N1*N2,DMIN,DMAX,COLBAR,NCOLS)
      TRCOL(6)=(DMAX-DMIN)/255.0
      TRCOL(4)=DMIN-TRCOL(6)
      IF (IPLOT.EQ.7) CALL LUTIN(LUTUSR,NCOLS,NCLUSR,IPLOT)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGVPORT(0.86,0.90,0.2,0.90)
      CALL PGWINDOW(1.0,2.0,DMIN,DMAX)
      CALL SETCOL(IPLOT,NCOLS,LUTUSR,RED,GR,BL)
      CALL PGSCH(0.50*SIZE)
      CALL PGBOX('B',0.0,0,'B',0.0,0)
      CALL PGBOX('C',0.0,0,'CMTSVI',0.0,0)
      CALL PGCELL(COLBAR,2,256,1,2,1,256,1.0,0.0,TRCOL,NCOLS,RED,GR,BL)
      CALL PGBOX('B',0.0,0,'B',0.0,0)
      CALL PGBOX('C',0.0,0,'CMTSVI',0.0,0)
      CALL PGSCH(SIZE)
      CALL PGVPORT(0.12,0.82,0.2,0.90)
      XMIN=TR(1)+TR(2)
      XMAX=TR(1)+FLOAT(NX)*TR(2)
      YMIN=TR(4)+TR(6)
      YMAX=TR(4)+FLOAT(NY)*TR(6)
      CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
      CALL PGBOX('C',0.0,0,'C',0.0,0)
      CALL PGBOX('BNSTI',0.0,0,'BNSTI',0.0,0)
      CALL PGCELL(MAP,N1,N2,1,NX,1,NY,1.0,0.0,TR,NCOLS,RED,GR,BL)
      CALL PGBOX('C',0.0,0,'C',0.0,0)
      CALL PGBOX('BNSTI',0.0,0,'BNSTI',0.0,0)
      END SUBROUTINE GREY

      SUBROUTINE STGREY(MAP,NMAP,FMIN,FMAX,COLBAR,NCOLS)
C     --------------------------------------------------
C
      REAL         MAP(*),COLBAR(2,*)
      CHARACTER*32 STRING
C
   1  WRITE(*,100) FMIN,FMAX
 100  FORMAT(' >>  Zmin & Zmax for plot ? (def=',1pe10.3,
     *       ',',E10.3,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) FMIN2,FMAX2
      IF (NNN.NE.0) THEN
        FMIN=FMIN2
        FMAX=FMAX2
      ENDIF
      IF (ABS(FMAX-FMIN).LT.1.0E-20) GOTO 1
      FNORM=1.0/(FMAX-FMIN)
   2  C=1.0
      WRITE(*,110) C
 110  FORMAT(' >>  Contrast factor ?  (def=',F3.1,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=2) C
      C=MAX(MIN(C,10.0),0.01)
      DO 10 I=1,NMAP
        F=(MAP(I)-FMIN)*FNORM
      IF (F.LE.0.0) THEN
          MAP(I)=0.0
        ELSEIF (F.GE.1.0) THEN
          MAP(I)=1.0
      ELSE
          MAP(I)=F**C
      ENDIF
  10  CONTINUE
      DCOL=0.999/FLOAT(NCOLS-1)
      COL=0.0
      IF (FMAX.LT.FMIN) THEN
        COL=1.0
        DCOL=-DCOL
      ENDIF
      DO 20 I=1,NCOLS
        COLBAR(1,I)=COL**C
        COLBAR(2,I)=COL**C
        COL=COL+DCOL
  20  CONTINUE
      END SUBROUTINE STGREY

      SUBROUTINE LUTIN(LUTUSR,NCOLS,NCLUSR,IPLOT)
C     -------------------------------------------
C
      CHARACTER*72 FILNAM
      REAL LUTUSR(3,NCOLS)
C
   1  WRITE(*,100)
 100  FORMAT(' INPUT> Filename for user colour-table ?  : ',$)
      READ(*,200,ERR=1) FILNAM
 200  FORMAT(A)
      OPEN(UNIT=17,FILE=FILNAM,STATUS='OLD',FORM='FORMATTED',ERR=1)
      DO 10 J=1,NCOLS
  10    READ(17,*,ERR=2,END=2) (LUTUSR(I,J),I=1,3)
   2  CLOSE(UNIT=17)
      NCLUSR=J-1
      WRITE(*,*)' No. of colour indicies read in = ',NCLUSR
      IF (NCLUSR.LE.1) THEN
        IPLOT=2
      ELSE
        NCOLS=NCLUSR
      ENDIF
      END SUBROUTINE LUTIN

      SUBROUTINE SETCOL(IPLOT,NCOLS,LUTUSR,R,G,B)
C     -------------------------------------------
C
      use coltabs
C      INCLUDE 'COLTABS.INC'
      REAL     LUTUSR(3,*),R(*),G(*),B(*)
C
      IF (IPLOT.EQ.2) THEN
        Z=0.0
        DZ=0.999/FLOAT(NCOLS-1)
        DO 10 I=1,NCOLS
          R(I)=Z
          G(I)=Z
          B(I)=Z
          Z=Z+DZ
  10    CONTINUE
      ELSEIF (IPLOT.EQ.3) THEN
        CALL STCOL1(HEAT,R,G,B,NCOLS)
      ELSEIF (IPLOT.EQ.4) THEN
        CALL STCOL1(SPECTRUM,R,G,B,NCOLS)
      ELSEIF (IPLOT.EQ.5) THEN
        CALL STCOL1(BGYRW,R,G,B,NCOLS)
      ELSEIF (IPLOT.EQ.6) THEN
        CALL STCOL1(SERP,R,G,B,NCOLS)
      ELSE
        CALL STCOL1(LUTUSR,R,G,B,NCOLS)
      ENDIF
      END SUBROUTINE SETCOL

      SUBROUTINE STCOL1(LUT,R,G,B,N)
C     ------------------------------
C
      REAL LUT(3,*),R(*),G(*),B(*)
C
      DO 10 I=1,N
        R(I)=LUT(1,I)
        G(I)=LUT(2,I)
        B(I)=LUT(3,I)
  10  CONTINUE
      END SUBROUTINE STCOL1

      SUBROUTINE VCOPY(X,Y,N)
C     -----------------------
C
      REAL X(*),Y(*)
C
      DO 10 I=1,N
  10    Y(I)=X(I)
      END SUBROUTINE VCOPY

      SUBROUTINE VRFILL(X,A,N)
C     ------------------------
C
      REAL X(*)
C
      DO 10 I=1,N
  10    X(I)=A
      END SUBROUTINE VRFILL

      SUBROUTINE LOGQYN(S,D,L)
C     ------------------------
C
      LOGICAL      L,L2
      CHARACTER*1  D,D2,A
      CHARACTER*45 STRING
      CHARACTER    S(*)
C
      IF (D.EQ.'Y') THEN
        L=.TRUE.
        D2='N'
        L2=.FALSE.
      ELSEIF (D.EQ.'N') THEN
        L=.FALSE.
        D2='Y'
        L2=.TRUE.
      ELSE
        WRITE(*,*)' Default should be Y or N !'
        RETURN
      ENDIF
      CALL STPACK(STRING,S,45)
   1  WRITE(*,100) STRING,D
 100  FORMAT(A45,' Y/N (def=',A1,')  : ',$)
      CALL FORMTQ(A,1,NNN)
      IF (NNN.EQ.0) RETURN
      IF (A.EQ.'y' .OR. A.EQ.'T' .OR. A.EQ.'t') A='Y'
      IF (A.EQ.'n' .OR. A.EQ.'F' .OR .A.EQ.'f') A='N'
      IF (A.EQ.D) THEN
        RETURN
      ELSEIF (A.EQ.D2) THEN
        L=L2
        RETURN
      ENDIF
      GOTO 1
      END SUBROUTINE LOGQYN

      SUBROUTINE STPACK(STRING,S,N)
C     -----------------------------
C
      CHARACTER STRING(*),S(*)
C
      DO 10 I=1,N
        STRING(I)=S(I)
        IF (S(I).EQ.'?') GOTO 20
  10  CONTINUE
  20  DO 30 J=I+1,N
  30    STRING(J)=' '
      END SUBROUTINE STPACK

      SUBROUTINE FORMTQ(STRING,NCHAR,NNN)
C     -----------------------------------
C
      CHARACTER*(*) STRING
C
      READ(*,200) STRING
 200  FORMAT(A)
      NNN=0
      DO 10 I=1,NCHAR
  10    IF (STRING(I:I).NE.' ') NNN=NNN+1
      END SUBROUTINE FORMTQ

      SUBROUTINE SBFINT(RGB,IC,IBMODE,IBUF,MAXBUF)
C     --------------------------------------------
C
      use grpckg1inc
      REAL             RGB(*)
C
      CHARACTER*16     TYPE,CHR
      LOGICAL          LPS,LCOLOR,LPSINT,PGNOTO
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
      REAL             RINIT(0:15),GINIT(0:15),BINIT(0:15)
      SAVE             LPSINT
      DATA             LPSINT /.FALSE./
      DATA     RINIT  /1.00,0.00,1.00,0.00,0.00,0.00,1.00,1.00,
     *                 1.00,0.50,0.00,0.00,0.50,1.00,0.33,0.67/
      DATA     GINIT  /1.00,0.00,0.00,1.00,0.00,1.00,0.00,1.00,
     *                 0.50,1.00,1.00,0.50,0.00,0.00,0.33,0.67/
      DATA     BINIT  /1.00,0.00,0.00,0.00,1.00,1.00,1.00,0.00,
     *                 0.00,0.00,0.50,1.00,1.00,0.50,0.33,0.67/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Initialises the software buffer for crystal-plotting. It should 
C    be called just once per plot (buffer), after PGWINDOW but before 
C    any crystal-related routines.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    RGB      R*4    I       3      The RGB values for the background.
C    IC       I*4    I       -      The index for the background colour.
C    IBMODE   I*4    I       -      Buffering mode for initialisation:
C                                     1 = Ordinary, default.
C                                     2 = Will want to save later.
C                                     3 = Initialise from saved buffers.
C    IBUF     I*4    I       -      Software buffer to be used (>=1).
C    MAXBUF   I*4    O       -      Maximum number of buffers available.
C
C Globals
C    SFTBUF
C    GRPCKG1.INC
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGNOTO     Logical function to test if a PGPLOT device is open.
C     GRBPIC     Sends a "begin picture" command to the device driver.
C     PGQINF     Inquires about general PGPLOT information.
C     GREXEC     Dispatches command to appropriate device driver.
C     COLINT     Sets up a colour table.
C     SBRFIL     Fills a real aray with a constant.
C     SBRCOP     Copies the contents of one real array to another.
C     SBFIN0     Inquires about viewport and window dimensions.
C
C History
C   D. S. Sivia       4 Apr 1995  Initial release.
C   D. S. Sivia      14 Sep 1995  Set up default colour indicies for PS.
C   D. S. Sivia      20 Oct 1995  Ignore NBUNCH and fix to one!
C   D. S. Sivia      15 Nov 1995  Allow initialisation to/from saved 
C                                 buffers.
C   D. S. Sivia       2 Aug 1996  Replaced pgplot.inc with SBFIN0, and
C                                 made appropriate additions to SFTBUF.
C   D. S. Sivia      16 Jul 1999  Added a couple of PGPLOT calls to 
C                                 force proper initialisation.
C-----------------------------------------------------------------------
C
C A PGPLOT initialisation precaution.
C
      IF (PGNOTO('SBFINT')) RETURN
      IF (.NOT.GRPLTD(GRCIDE)) CALL GRBPIC
C
C Begin SBFINT proper.
C
      CALL SBFIN0(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *            YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
      LPS=.TRUE.
      NDOTS=100
      CALL PGQINF('TYPE',TYPE,LCHR)
      IF (TYPE.EQ.'PS' .OR. TYPE.EQ.'VPS') THEN
        LCOLOR=.FALSE.
        NXP=NINT(FLOAT(NDOTS)*XLEN/XPERIN)
        NYP=NINT(FLOAT(NDOTS)*YLEN/YPERIN)
      ELSEIF (TYPE.EQ.'CPS' .OR. TYPE.EQ.'VCPS') THEN
        LCOLOR=.TRUE.
        NXP=NINT(FLOAT(NDOTS)*XLEN/XPERIN)
        NYP=NINT(FLOAT(NDOTS)*YLEN/YPERIN)
      ELSE
        LPS=.FALSE.
        NBUF=0
        LCHR=LEN(CHR)
        CALL GREXEC(GRGTYP,4,RBUF,NBUF,CHR,LCHR)
        IF (CHR(7:7).EQ.'P') THEN
          NXP=1+NINT(XLEN)
          NYP=1+NINT(YLEN)
        ELSE
          STOP' Sorry, SFBINT does not support this device !'
        ENDIF
      ENDIF
      IBFMOD=IBMODE
      NTOT=NXP*NYP
      MAXBUF=INT(2000000/NTOT)-1
      IF (IBFMOD.EQ.2 .OR. IBFMOD.EQ.3) MAXBUF=MAXBUF-2 
      KSTART=NTOT*MAX(MIN(IBUF,MAXBUF),1)
      IF ((KSTART+NTOT).LE.2000000 .AND. MAXBUF.GE.1) THEN
        IF (IBFMOD.EQ.3) THEN
          IZSAVE=NTOT*(1+MAXBUF)+1
          ICSAVE=IZSAVE+NTOT
          CALL SBRCOP(SBBUFF(IZSAVE),SBBUFF(1),NTOT)
          CALL SBRCOP(SBBUFF(ICSAVE),SBBUFF(KSTART+1),NTOT)
        ELSE
          CALL SBRFIL(SBBUFF(1),-1.0E20,NTOT)
          CALL SBRFIL(SBBUFF(KSTART+1),FLOAT(IC),NTOT)
        ENDIF
      ELSE
        STOP' Sorry, the software buffer is too small !'
      ENDIF
      CALL COLINT(RGB,IC,IC,0.5,0.0,1.0)
      IF (LPS .AND. .NOT. LPSINT) THEN
        DO 10 I=0,15
          IF (LCOLOR) THEN
            IR(I)=NINT(255.0*RINIT(I))
            IG(I)=NINT(255.0*GINIT(I))
            IB(I)=NINT(255.0*BINIT(I))
          ELSE
            IR(I)=NINT(85.0*(RINIT(I)+GINIT(I)+BINIT(I)))
          ENDIF
  10    CONTINUE
        LPSINT=.TRUE.
      ENDIF
      END SUBROUTINE SBFINT

      SUBROUTINE SBFIN0(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                  YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
C     ---------------------------------------------------------
C
      DATA BIG,SMALL /1.0E+20,1.0E-20/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Obtain some viewport and window information about the current 
C    PGPLOT device, without directly accessing the common blocks in
C    pgplot.inc.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    XPERIN   R*4    O       -      Plot X scale in dots/inch.
C    YPERIN   R*4    O       -      Plot Y scale in dots/inch.
C    XOFF     R*4    O       -      Absolute coord of blc of viewport.
C    YOFF     R*4    O       -      Absolute coord of blc of viewport.
C    XLEN     R*4    O       -      Width of viewport in absolute coord.
C    YLEN     R*4    O       -      Height of viewport in absolute coord.
C    XORG     R*4    O       -      Absolute coord of world X=0.
C    YORG     R*4    O       -      Absolute coord of world Y=0.
C    XSCALE   R*4    O       -      Absolute units per world coord in X.
C    YSCALE   R*4    O       -      Absolute units per world coord in Y.
C    XBLC     R*4    O       -      World X coord at blc of window.
C    XTRC     R*4    O       -      World X coord at trc of window.
C    YBLC     R*4    O       -      World Y coord at blc of window.
C    YTRC     R*4    O       -      World Y coord at trc of window.
C
C Globals
C     None.
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGQVP      Inquires about viewport dimensions.
C     PGQWIN     Inquires about world coords of window.
C
C History
C   D. S. Sivia       2 Aug 1996  Initial release.
C-----------------------------------------------------------------------
C
      CALL PGQWIN(XBLC,XTRC,YBLC,YTRC)
      CALL PGQVP(1,XI1,XI2,YI1,YI2)
      CALL PGQVP(3,XOFF,XP2,YOFF,YP2)
      XLEN=ABS(XP2-XOFF)
      YLEN=ABS(YP2-YOFF)
      XPERIN=XLEN/(ABS(XI2-XI1)+SMALL)
      YPERIN=YLEN/(ABS(YI2-YI1)+SMALL)
      XWDIF=XTRC-XBLC
      YWDIF=YTRC-YBLC
      AXWDIF=BIG
      AYWDIF=BIG
      IF (ABS(XWDIF).GT.SMALL) AXWDIF=1.0/XWDIF
      IF (ABS(YWDIF).GT.SMALL) AYWDIF=1.0/YWDIF
      XSCALE=XLEN*AXWDIF
      YSCALE=YLEN*AYWDIF
      XORG=(XOFF*XTRC-XP2*XBLC)*AXWDIF
      YORG=(YOFF*YTRC-YP2*YBLC)*AYWDIF
      END SUBROUTINE SBFIN0

      SUBROUTINE SBFSAV(IBUF)
C     -----------------------
C
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Save a rendered picture-buffer, and its Z-buffer, for subsequent 
C    use in re-initialisation with SBFINT.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    IBUF     I*4    I       -      Software buffer to be saved (>=1).
C
C Globals
C    SFTBUF
C    GRPCKG1.INC
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBRCOP     Copies the contents of one real array to another.
C
C History
C   D. S. Sivia      15 Nov 1995  Initial release.
C-----------------------------------------------------------------------
C
      NTOT=NXP*NYP
      MAXBUF=INT(2000000/NTOT)-3
      IBUF1=MAX(IBUF,1)
      IF (IBUF1.GT.MAXBUF) RETURN
      KSTART=NTOT*IBUF1
      IZSAVE=NTOT*(1+MAXBUF)+1
      ICSAVE=IZSAVE+NTOT
      CALL SBRCOP(SBBUFF(1),SBBUFF(IZSAVE),NTOT)
      CALL SBRCOP(SBBUFF(KSTART+1),SBBUFF(ICSAVE),NTOT)
      END SUBROUTINE SBFSAV

      SUBROUTINE COLINT(RGB,IC1,IC2,DIFUSE,SHINE,POLISH)
C     --------------------------------------------------
C
      REAL             RGB(*)
C
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Initialises a colour table for a geometrical object. In general,
C    it is recommended that SHINE = 0.0 if DIFUSE > 0.0 and vice versa.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    RGB      R*4    I       3      Red, green and blue intenisty for 
C                                   fully-lit non-shiny object (0-1).
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for shading.
C    DIFUSE   R*4    I       -      Diffusiveness of object (0-1).
C    SHINE    R*4    I       -      Whiteness of bright spot (0-1).
C    POLISH   R*4    I       -      Controls size of bright spot.
C
C Globals
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGQCOL     Inquires about the colour capability.
C     PGSCR      Assigns an RGB value to a colour index.
C
C History
C   D. S. Sivia       4 Apr 1995  Initial release.
C-----------------------------------------------------------------------
C
      IF (RGB(1).LT.0.0 .OR. RGB(1).GT.1.0) RETURN
      IF (RGB(2).LT.0.0 .OR. RGB(2).GT.1.0) RETURN
      IF (RGB(3).LT.0.0 .OR. RGB(3).GT.1.0) RETURN
      IF (DIFUSE.LT.0.0 .OR. DIFUSE.GT.1.0) RETURN
      IF (SHINE.LT.0.0 .OR. SHINE.GT.1.0) RETURN
      IF (IC2.LT.IC1) RETURN
      IF (LPS) THEN
        ICMIN=0
        ICMAX=255
      ELSE
        CALL PGQCOL(ICMIN,ICMAX)
      ENDIF
      IF (IC1.LT.ICMIN .OR. IC2.GT.ICMAX) THEN
        WRITE(*,*)' Invalid colour-indicies for the chosen device !'
        RETURN
      ENDIF
      POLSH2=MAX(POLISH/2.0,0.5)
      NC=IC2-IC1
      IF (NC.EQ.0) THEN
        RED=RGB(1)
        GRN=RGB(2)
        BLU=RGB(3)
      ELSE
        GREY=0.0
        DGREY=1.0/FLOAT(NC)
        DR=DIFUSE*RGB(1)/FLOAT(NC)
        DG=DIFUSE*RGB(2)/FLOAT(NC)
        DB=DIFUSE*RGB(3)/FLOAT(NC)
        RED=MAX(RGB(1)*(1.0-DIFUSE),0.0)
        GRN=MAX(RGB(2)*(1.0-DIFUSE),0.0)
        BLU=MAX(RGB(3)*(1.0-DIFUSE),0.0)
        R=RED
        G=GRN
        B=BLU
      ENDIF
      DO 10 IC=IC1,IC2
        IF (.NOT. LPS) THEN
          CALL PGSCR(IC,RED,GRN,BLU)
        ELSE
          IF (LCOLOR) THEN
            IR(IC)=NINT(255.0*RED)
            IG(IC)=NINT(255.0*GRN)
            IB(IC)=NINT(255.0*BLU)
          ELSE
            IR(IC)=NINT(85.0*(RED+GRN+BLU))
          ENDIF
        ENDIF
        R=R+DR
        G=G+DG
        B=B+DB
        GREY=GREY+DGREY
        POLSHN=SHINE*(GREY**POLSH2)
        RED=MIN(R+POLSHN,1.0)
        GRN=MIN(G+POLSHN,1.0)
        BLU=MIN(B+POLSHN,1.0)
  10  CONTINUE
      END SUBROUTINE COLINT

      SUBROUTINE COLTAB(RGB,NCOL,ALFA,IC1,IC2)
C     ----------------------------------------
C
      REAL             RGB(3,*)
C
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Initialises a colour table for a "grey-scale" map.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    RGB      R*4    I   3 X NCOL   Red, green and blue intenisty for 
C                                   the colour table.
C    NCOL     I*4    I       -      No. of colours in the input table.
C    ALFA     R*4    I       -      Contrast-factor (linear=1).
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the output.
C
C Globals
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGQCOL     Inquires about the colour capability.
C     PGSCR      Assigns an RGB value to a colour index.
C
C History
C   D. S. Sivia      30 Apr 1995  Initial release.
C-----------------------------------------------------------------------
C
      IF (IC2.LE.IC1) RETURN
      IF (LPS) THEN
        ICMIN=0
        ICMAX=255
      ELSE
        CALL PGQCOL(ICMIN,ICMAX)
      ENDIF
      IF (IC1.LT.ICMIN .OR. IC2.GT.ICMAX) THEN
        WRITE(*,*)' Invalid colour-indicies for the chosen device !'
        RETURN
      ENDIF
      NC=IC2-IC1
      COL=0.0
      DCOL=1.0/FLOAT(NC)
      DO 10 I=0,NC
        COLALF=FLOAT(NCOL-1)*MIN(COL**ALFA,0.99999)
        ICOL=INT(COLALF)
        DICOL=COLALF-FLOAT(ICOL)
        RL=RGB(1,ICOL+1)+DICOL*(RGB(1,ICOL+2)-RGB(1,ICOL+1))
        GL=RGB(2,ICOL+1)+DICOL*(RGB(2,ICOL+2)-RGB(2,ICOL+1))
        BL=RGB(3,ICOL+1)+DICOL*(RGB(3,ICOL+2)-RGB(3,ICOL+1))
        R=MIN(MAX(RL,0.0),1.0)
        G=MIN(MAX(GL,0.0),1.0)
        B=MIN(MAX(BL,0.0),1.0)
        IF (LPS) THEN
          IF (LCOLOR) THEN
            IR(IC1+I)=NINT(R*255.0)
            IG(IC1+I)=NINT(G*255.0)
            IB(IC1+I)=NINT(B*255.0)
          ELSE
            IR(IC1+I)=NINT((R+G+B)*85.0)
          ENDIF
        ELSE
          CALL PGSCR(IC1+I,R,G,B)
        ENDIF
        COL=COL+DCOL
  10  CONTINUE
      END SUBROUTINE COLTAB

      SUBROUTINE COLSRF(RGB,NCOL,ALFA,IC1,IC2,NCBAND,DIFUSE,SHINE,
     *                  POLISH)
C     ------------------------------------------------------------
C
      REAL             RGB(3,*)
C
      REAL             RGBOUT(3)
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Initialises a colour table for a 3-D surface rendering of a 2-D 
C   array of "data".
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    RGB      R*4    I   3 X NCOL   Red, green and blue intenisty for 
C                                   the colour table.
C    NCOL     I*4    I       -      No. of colours in the input table.
C    ALFA     R*4    I       -      Contrast-factor (linear=1).
C    IC1,IC2  I*4    I       -      Lowest and highest colour-index to
C                                   be used for the rendering.
C    NCBAND   I*4    I       -      Number of colour-bands for the
C                                   height, so that the number of shades
C                                   per band = (IC2-IC1+1)/NCBAND.
C    DIFUSE   R*4    I       -      Diffusiveness of object (0-1).
C    SHINE    R*4    I       -      Whiteness of bright spot (0-1).
C    POLISH   R*4    I       -      Controls size of bright spot.
C
C Globals
C    SFTBUF
C
C External Calls
!




C   SUBROUTINE   DESCRIPTION
C     PGQCOL     Inquires about the colour capability.
C
C History
C   D. S. Sivia      30 Oct 1995  Initial release.
C-----------------------------------------------------------------------
C
      IF (IC2.LE.IC1) RETURN
      IF (LPS) THEN
        ICMIN=0
        ICMAX=255
      ELSE
        CALL PGQCOL(ICMIN,ICMAX)
      ENDIF
      IF (IC1.LT.ICMIN .OR. IC2.GT.ICMAX) THEN
        WRITE(*,*)' Invalid colour-indicies for the chosen device !'
        RETURN
      ENDIF
      NSHADS=MAX((IC2-IC1+1)/MAX(NCBAND,1),1)
      COL=0.0
      DCOL=1.0/FLOAT(MAX(NCBAND-1,1))
      IC=IC1
      DO 10 I=1,NCBAND
        COLALF=FLOAT(NCOL-1)*MIN(COL**ALFA,0.99999)
        ICOL=INT(COLALF)
        DICOL=COLALF-FLOAT(ICOL)
        RL=RGB(1,ICOL+1)+DICOL*(RGB(1,ICOL+2)-RGB(1,ICOL+1))
        GL=RGB(2,ICOL+1)+DICOL*(RGB(2,ICOL+2)-RGB(2,ICOL+1))
        BL=RGB(3,ICOL+1)+DICOL*(RGB(3,ICOL+2)-RGB(3,ICOL+1))
        RGBOUT(1)=MIN(MAX(RL,0.0),1.0)
        RGBOUT(2)=MIN(MAX(GL,0.0),1.0)
        RGBOUT(3)=MIN(MAX(BL,0.0),1.0)
        CALL COLINT(RGBOUT,IC,IC+NSHADS-1,DIFUSE,SHINE,POLISH)
        IC=IC+NSHADS
        COL=COL+DCOL
  10  CONTINUE
      END SUBROUTINE COLSRF

      SUBROUTINE SBFCLS(IBUF)
C     -----------------------
C
      use grpckg1inc
      REAL            BUFFER(2050),RNDBUF(0:3050)
      INTEGER         VALUE(33)
      LOGICAL         LPS,LCOLOR,LSTART
      CHARACTER       INLINE*80,CHR*16
      SAVE            LSTART,RNDBUF
      COMMON /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
      DATA            LODD,RNDSCL,RNDOFF,NRND /37614625,0.96,-0.5,3000/
      DATA            LSTART /.FALSE./
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Closes the software buffer for crystal-plotting, by outputing it
C    to the screen or writing out a postscript file (as appropriate).
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    IBUF     I*4    I       -      Software buffer to be output (>=1).
C
C Globals
C    SFTBUF
C    GRPCKG1.INC
C
C External Calls
!




C   SUBROUTINE   DESCRIPTION
C     GREXEC     Dispatches command to appropriate device driver.
C     PGBBUF     Recommended initial call (to start a PGPLOT buffer).
C     PGEBUF     Recommended final call (to end a PGPLOT buffer).
C
C History
C   D. S. Sivia       4 Apr 1995  Initial release.
C   D. S. Sivia      20 Oct 1995  Hardwire assumption that NGANG=1.
C   D. S. Sivia      21 Oct 1997  Made slightly friendlier for NT.
C   D. S. Sivia      30 Jan 1998  Explicitly SAVE array RNDBUF.
C-----------------------------------------------------------------------
C
      NTOT=NXP*NYP
      MAXBUF=INT(2000000/NTOT)-1
      IF (IBFMOD.EQ.2 .OR. IBFMOD.EQ.3) MAXBUF=MAXBUF-2
      IF (IBUF.GT.MAXBUF) RETURN
      KEND=NTOT*MAX(IBUF,1)
      IF (.NOT. LSTART) THEN
        DO 10 I=0,NRND-1
#ifdef CRAY
          CALL RANDOM_NUMBER(HARVEST=RAND)
  10      RNDBUF(I)=RNDSCL*(RAND+RNDOFF)
#elif RISC
  10      RNDBUF(I)=RNDSCL*(erand48(LODD)+RNDOFF)
#else
  10      RNDBUF(I)=RNDSCL*(RAN(LODD)+RNDOFF)
#endif
        LSTART=.TRUE.
      ENDIF
      IF (LPS) THEN
        AT=FLOAT(NXP)/XLEN
        BT=0.0
        CT=0.0
        DT=FLOAT(NYP)/YLEN
        TX=-AT*XOFF
        TY=-DT*YOFF
        WRITE (INLINE,100) NINT(XOFF),NINT(YOFF),NINT(XLEN),NINT(YLEN),
     *                    -NINT(XLEN)
 100    FORMAT(I6,I6,' moveto ',I6, ' 0 rlineto  0 ',I6,' rlineto ',
     *         I6,' 0 rlineto')
        CALL GRTERM
        CALL GRESC(' newpath ')
        CALL GRESC(INLINE)
        CALL GRESC(' closepath ')
        WRITE (INLINE, '(A,I5,A)') '/picstr ',NXP,' string def'
        CALL GRESC(INLINE)
        WRITE (INLINE,110) NXP,NYP,AT,BT,CT,DT,TX,TY
 110    FORMAT(2I4,' 8 [',6(1PE10.3,' '),']')
        CALL GRESC(INLINE)
        CALL GRESC('{ currentfile picstr readhexstring pop}')
        IF (LCOLOR) THEN
          CALL GRESC('  false 3 colorimage')
        ELSE
          CALL GRESC('  image')
        ENDIF
        L=0
        DO 20 K=KEND+1,KEND+NTOT
          L=L+1
          IC=NINT(SBBUFF(K)+RNDBUF(MOD(K,NRND)))
          IF (LCOLOR) THEN
            VALUE(L)=IR(IC)
            VALUE(L+1)=IG(IC)
            VALUE(L+2)=IB(IC)
            L=L+2
          ELSE
            VALUE(L)=IR(IC)
          ENDIF
          IF (L.EQ.33) THEN
            WRITE(INLINE,120) (VALUE(I),I=1,33)
 120        FORMAT(33Z2.2)
            CALL GRESC(INLINE(1:66))
            L=0
          ENDIF
  20    CONTINUE
        IF (L.NE.0) THEN
          WRITE(INLINE,120) (VALUE(I),I=1,L)
          CALL GRESC(INLINE(1:2*L))
        ENDIF
        CALL GRESC(' newpath ')
        CALL GRTERM
      ELSE
        CALL PGBBUF
        NXP2=NXP+2
        K=KEND+1
        RND900=900.0
        BUFFER(1)=FLOAT(NINT(XOFF))
        BUFFER(2)=FLOAT(NINT(YOFF))
        DO 40 J=1,NYP
          L=1+INT(RND900*(RNDBUF(J)-RNDOFF))
          CALL SBCLS0(BUFFER(3),SBBUFF(K),RNDBUF(L),NXP)
          CALL GREXEC(GRGTYP,26,BUFFER,NXP2,CHR,LCHR)
          K=K+NXP
          BUFFER(2)=BUFFER(2)+1.0
  40    CONTINUE
        CALL PGEBUF
      ENDIF
      END SUBROUTINE SBFCLS

      SUBROUTINE SBCLS0(A,B,C,N)
C     --------------------------
C
      REAL A(*),B(*),C(*)
C
      DO 10 I=1,N
  10    A(I)=FLOAT(NINT(B(I)+C(I)))
      END SUBROUTINE SBCLS0

      SUBROUTINE SBBALL(EYE,CENTRE,RADIUS,IC1,IC2,LIGHT,LSHINE,X0,Y0,R0)
C     ------------------------------------------------------------------
C
      REAL            EYE(*),CENTRE(*),LIGHT(*)
      LOGICAL         LSHINE

      REAL            SURF(3)
      REAL*8          ALFA,BETA,GAMA,XMU,A,B,C,DET,Q,DX0H,DY0H
      REAL*8          DZE,DZE2,DGAMZE,DBL1,DBL2,DSMALL
      REAL*8          XL0,XL1,HYP,SINPHI,COSPHI,R1,R2
      LOGICAL         LPS,LCOLOR
      COMMON /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots a shiny or matt coloured ball. All 
C    (x,y,z) values are taken to be given in world coordinates. The 
C    z-component of the eye-poisition should be positive and that of
C    the ball-centre should be negative (< -radius); the viewing-screen
C    is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    CENTRE   R*4    I       3      (x,y,z) coordinate of ball-centre.
C    RADIUS   R*4    I       -      Radius of ball.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    LSHINE   L*1    I       -      Shiny ball if .TRUE., else diffuse.
C    X0,Y0    R*4    O       -      Centre of projected ball.
C    R0       R*4    O       -      Average radius of projected ball.
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBGLOS     Works out colour-shade for surface of ball.
C
C History
C   D. S. Sivia       7 Apr 1995  Initial release.
C-----------------------------------------------------------------------
C
C Some initial checks.
C
      SMALL=1.0E-20
      DSMALL=DBLE(SMALL)
      IF (EYE(3).LE.0.0) RETURN
      IF (RADIUS.LE.0.0) RETURN
      IF (CENTRE(3).GT.-RADIUS) RETURN
C
C Calculate parameters of projected ellipse.
C
      DZE=DBLE(EYE(3))
      DZE2=DZE**2
      ALFA=DBLE(EYE(1)-CENTRE(1))
      BETA=DBLE(EYE(2)-CENTRE(2))
      GAMA=DBLE(EYE(3)-CENTRE(3))
      XMU=DBLE(RADIUS**2)-(ALFA**2+BETA**2+GAMA**2)
      A=XMU*(XMU+ALFA**2)
      B=XMU*ALFA*BETA
      C=XMU*(XMU+BETA**2)
      DET=ABS(A*C-B**2)+DSMALL
      DX0H=GAMA*XMU*DZE*(ALFA*C-BETA*B)/DET
      DY0H=GAMA*XMU*DZE*(BETA*A-ALFA*B)/DET
      Q=A*DX0H**2+2.0*B*DX0H*DY0H+C*DY0H**2-XMU*(XMU+GAMA**2)*DZE2
      X0H=SNGL(DX0H)
      Y0H=SNGL(DY0H)
      DX=(XTRC-XBLC)/FLOAT(NXP-1)
      DY=(YTRC-YBLC)/FLOAT(NYP-1)
      XDIF=SNGL(SQRT(ABS(C*Q/DET)+DSMALL))
      XMIN=X0H-XDIF+EYE(1)
      XMAX=X0H+XDIF+EYE(1)
      IXMIN=INT((XMIN-XBLC)/DX)+2
      IXMAX=INT((XMAX-XBLC)/DX)+1
      IF (IXMIN.GT.NXP .OR. IXMAX.LT.1) RETURN
      YDIF=(SQRT(ABS(A*Q/DET)+DSMALL))
      YMIN=Y0H-YDIF+EYE(2)
      YMAX=Y0H+YDIF+EYE(2)
      JYMIN=INT((YMIN-YBLC)/DY)+2
      JYMAX=INT((YMAX-YBLC)/DY)+1
      IF (JYMIN.GT.NYP .OR. JYMAX.LT.1) RETURN
      IF (JYMIN.LT.1) JYMIN=1
      IF (JYMAX.GT.NYP) JYMAX=NYP
      ZMAX=CENTRE(3)+RADIUS
      X0=X0H+EYE(1)
      Y0=Y0H+EYE(2)
      COREL=SNGL(SQRT(ABS((B*B)/(A*C))+DSMALL))
      IF (COREL.GT.0.0001) THEN
        XL0=(A+C)/2.0D0
        XL1=XL0-SQRT(ABS(XL0*XL0-DET)+DSMALL)
        HYP=SQRT((XL1-A)**2+B**2+DSMALL)
        SINPHI=(XL1-A)/HYP
        COSPHI=B/HYP
      ELSE
        SINPHI=0.0D0
        COSPHI=1.0D0
      ENDIF
      R1=SQRT(Q/(A*COSPHI*COSPHI+SINPHI*(C*SINPHI+2.0*B*COSPHI)))
      R2=SQRT(Q/(A*SINPHI*SINPHI+COSPHI*(C*COSPHI-2.0*B*SINPHI)))
      R0=SNGL((R1+R2)/2.0D0)
C
C Fill the inside of the projected ellipse with the right colours.
C
      NC=IC2-IC1
      COL0=FLOAT(IC1)
      COLSCL=FLOAT(NC)
      XL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
      XN2=RADIUS**2
      XNL2=1.0/SQRT(XN2*XL2+SMALL)
      XN3=1.0/(XN2+SMALL)
      YH=YBLC+DY*FLOAT(JYMIN-1)-EYE(2)
      DGAMZE=GAMA*DZE
      BA=SNGL(B/A)
      DO 20 JY=JYMIN,JYMAX
        YH2=YH**2
        BETAYH=BETA*YH
        YDIF=YH-Y0H
        XDIF=SNGL(SQRT(ABS(A*Q-DET*DBLE(YDIF**2))+DSMALL)/A)
        XMIN=X0H-BA*YDIF-XDIF+EYE(1)
        XMAX=X0H-BA*YDIF+XDIF+EYE(1)
        IXMIN=INT((XMIN-XBLC)/DX)+2
        IXMAX=INT((XMAX-XBLC)/DX)+1
        IF (IXMIN.LE.NXP .AND. IXMAX.GE.1) THEN
          IF (IXMIN.LT.1) IXMIN=1
          IF (IXMAX.GT.NXP) IXMAX=NXP
          XH=XBLC+DX*FLOAT(IXMIN-1)-EYE(1)
          K=(JY-1)*NXP+IXMIN
          DO 10 IX=IXMIN,IXMAX
            IF (ZMAX.GT.SBBUFF(K)) THEN
              XH2=XH**2
              ALFAXH=ALFA*XH
              DBL1=DBLE(ALFAXH+BETAYH)-DGAMZE
              DBL2=DBLE(XH2+YH2)+DZE2
              XLM=SNGL((-DBL1-SQRT(ABS(DBL1**2+XMU*DBL2)+DSMALL))/DBL2)
              SURF(3)=EYE(3)*(1.0-XLM)
              IF (SURF(3).GT.SBBUFF(K)) THEN
                SBBUFF(K)=SURF(3)
                IF (NC.EQ.0) THEN
                  SBBUFF(KSTART+K)=COL0
                ELSE
                  SURF(2)=EYE(2)+YH*XLM
                  SURF(1)=EYE(1)+XH*XLM
                  CALL SBGLOS(EYE,CENTRE,LIGHT,SURF,XNL2,XN3,SMALL,
     *                        LSHINE,COLOUR)
                  SBBUFF(KSTART+K)=COL0+COLOUR*COLSCL
                ENDIF
              ENDIF
            ENDIF
            K=K+1
            XH=XH+DX
  10      CONTINUE
        ENDIF
        YH=YH+DY
  20  CONTINUE
      END SUBROUTINE SBBALL

      SUBROUTINE SBGLOS(EYE,CENTRE,LIGHT,SURF,XNL2,XN3,SMALL,LSHINE,
     *                  COLOUR)
C     --------------------------------------------------------------
C
!




C Support subroutine for SBBALL, to work out colour-shade.
C
      REAL    EYE(*),CENTRE(*),LIGHT(*),SURF(*)
      LOGICAL LSHINE
      REAL    NORMAL(3),REFLEC(3),VIEW(3)
C
      COLOUR=0.0
      XNL=0.0
      DO 10 I=1,3
        NORMAL(I)=SURF(I)-CENTRE(I)
        XNL=XNL+NORMAL(I)*LIGHT(I)
  10  CONTINUE
      IF (XNL.GE.0.0) RETURN
      IF (LSHINE) THEN
        RFNORM=(XNL+XNL)*XN3
        XRV=0.0
        DO 20 I=1,3
          VIEW(I)=EYE(I)-SURF(I)
          REFLEC(I)=LIGHT(I)-RFNORM*NORMAL(I)
          XRV=XRV+REFLEC(I)*VIEW(I)
  20    CONTINUE
        IF (XRV.LT.0.0) RETURN
        REF2=0.0
        VIEW2=0.0
        DO 30 I=1,3
          REF2=REF2+REFLEC(I)**2
          VIEW2=VIEW2+VIEW(I)**2
  30    CONTINUE
        COLOUR=MIN(XRV**2/(ABS(REF2*VIEW2)+SMALL),1.0)
      ELSE
        COLOUR=MIN(-XNL*XNL2,1.0)
      ENDIF
      END SUBROUTINE SBGLOS

      SUBROUTINE SBLINE(EYE,END1,END2,ICOL,LDASH)
C     -------------------------------------------
C
      REAL             EYE(*),END1(*),END2(*)
      LOGICAL          LDASH
C
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
      DATA             NDASH1,NDASH2 /7,3/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
!




C      This subroutine draws a straight line between two points. All 
C    (x,y,z) values are taken to be given in world coordinates. The 
C    z-component of the eye-poisition should be positive, while that 
C    of both the ends should be negative; the viewing-screen is fixed
C    at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    END1     R*4    I       3      (x,y,z) coordinate of end-1.
C    END2     R*4    I       3      (x,y,z) coordinate of end-2.
C    ICOL     I*4    I       -      Colour-index for line.
C    LDASH    L*1    I       -      Dashed line if .TRUE. (else cont.).
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBLIN1     Calculates the projection of (x,y,z) on viewing screen.
C
C History
C   D. S. Sivia       4 Apr 1995  Initial release.
C   D. S. Sivia      25 Oct 1997  Prevent occasional Z-coordinate bug.
C-----------------------------------------------------------------------
C
C Some initial checks and clipping.
C
      SMALL=1.0E-10
      IF (EYE(3).LE.0.0) RETURN
      IF (END1(3).GE.0.0 .OR. END2(3).GE.0.0) RETURN
      CALL SBLIN1(EYE,END1(1),END1(2),END1(3),XW1,YW1)
      CALL SBLIN1(EYE,END2(1),END2(2),END2(3),XW2,YW2)
      XDIF=XW2-XW1
      YDIF=YW2-YW1
      IF ((ABS(XDIF)+ABS(YDIF)).LT.SMALL) RETURN
      IF (ABS(XDIF).LT.SMALL) XDIF=SMALL
      IF (ABS(YDIF).LT.SMALL) YDIF=SMALL
      IF (XW1.LE.XBLC) THEN
        IF (XW2.LE.XBLC) RETURN
        XLAM=(XBLC-XW1)/XDIF
        XW1=XW1+XLAM*XDIF
        YW1=YW1+XLAM*YDIF
      ELSEIF (XW1.GE.XTRC) THEN
        IF (XW2.GE.XTRC) RETURN
        XLAM=(XTRC-XW1)/XDIF
        XW1=XW1+XLAM*XDIF
        YW1=YW1+XLAM*YDIF
      ENDIF
      IF (YW1.LE.YBLC) THEN
        IF (YW2.LE.YBLC) RETURN
        YLAM=(YBLC-YW1)/YDIF
        XW1=XW1+YLAM*XDIF
        YW1=YW1+YLAM*YDIF
      ELSEIF (YW1.GE.YTRC) THEN
        IF (YW2.GE.YTRC) RETURN
        YLAM=(YTRC-YW1)/YDIF
        XW1=XW1+YLAM*XDIF
        YW1=YW1+YLAM*YDIF
      ENDIF
      IF (XW2.LE.XBLC) THEN
        XLAM=(XBLC-XW1)/XDIF
        XW2=XW1+XLAM*XDIF
        YW2=YW1+XLAM*YDIF
      ELSEIF (XW2.GE.XTRC) THEN
        XLAM=(XTRC-XW1)/XDIF
        XW2=XW1+XLAM*XDIF
        YW2=YW1+XLAM*YDIF
      ENDIF
      IF (YW2.LE.YBLC) THEN
        YLAM=(YBLC-YW1)/YDIF
        XW2=XW1+YLAM*XDIF
        YW2=YW1+YLAM*YDIF
      ELSEIF (YW2.GE.YTRC) THEN
        YLAM=(YTRC-YW1)/YDIF
        XW2=XW1+YLAM*XDIF
        YW2=YW1+YLAM*YDIF
      ENDIF
C
      COLOUR=FLOAT(ICOL)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      X1=1.0+(XW1-XBLC)*DX
      X2=1.0+(XW2-XBLC)*DX
      XDIF=X2-X1
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      Y1=1.0+(YW1-YBLC)*DY
      Y2=1.0+(YW2-YBLC)*DY
      YDIF=Y2-Y1
      ZW1=END1(3)
      IF (ABS(XDIF).GE.ABS(YDIF)) THEN
        IF (ABS(XDIF).LT.SMALL) RETURN
        IF (X1.LE.X2) THEN
          A=EYE(3)-END1(3)
          B=EYE(3)*(EYE(1)-END1(1))
          C=END2(3)-END1(3)
          D=EYE(3)*(END2(1)-END1(1))
        ELSE
          X3=X1
          Y3=Y1
          X1=X2
          Y1=Y2
          X2=X3
          Y2=Y3
          ZW1=END2(3)
          A=EYE(3)-END2(3)
          B=EYE(3)*(EYE(1)-END2(1))
          C=END1(3)-END2(3)
          D=EYE(3)*(END1(1)-END2(1))
        ENDIF
        IF (ABS(C).LT.SMALL) C=SMALL
        D=D/C
        IF (ABS(D).LT.SMALL) D=SMALL
        DYJ=YDIF/XDIF
        IX1=NINT(X1)
        IX2=NINT(X2)
        YJ=Y1+DYJ*(FLOAT(IX1)-X1)
        DXX=1.0/DX
        IF (LDASH) THEN
          DYJJ=(NDASH1-NDASH2-1)*DYJ
          DO 20 II=IX1,IX2,NDASH1
            III=MIN(II+NDASH2,IX2)
            XX=XBLC-EYE(1)+DXX*FLOAT(II-1)
            DO 10 I=II,III
              K=NXP*(NINT(YJ)-1)+I
              Z=ZW1+(A*XX+B)/(XX+D)
              IF (Z.GT.SBBUFF(K)) THEN
                SBBUFF(K)=Z
                SBBUFF(KSTART+K)=COLOUR
              ENDIF
              YJ=YJ+DYJ
              XX=XX+DXX
  10        CONTINUE
            YJ=YJ+DYJJ
  20      CONTINUE
        ELSE
          XX=XBLC-EYE(1)+DXX*FLOAT(IX1-1)
          DO 30 I=IX1,IX2
            K=NXP*(NINT(YJ)-1)+I
            Z=ZW1+(A*XX+B)/(XX+D)
            IF (Z.GT.SBBUFF(K)) THEN
              SBBUFF(K)=Z
              SBBUFF(KSTART+K)=COLOUR
            ENDIF
            YJ=YJ+DYJ
            XX=XX+DXX
  30      CONTINUE
        ENDIF
      ELSE
        IF (ABS(YDIF).LT.SMALL) RETURN
        IF (Y1.LE.Y2) THEN
          A=EYE(3)-END1(3)
          B=EYE(3)*(EYE(2)-END1(2))
          C=END2(3)-END1(3)
          D=EYE(3)*(END2(2)-END1(2))
        ELSE
          X3=X1
          Y3=Y1
          X1=X2
          Y1=Y2
          X2=X3
          Y2=Y3
          ZW1=END2(3)
          A=EYE(3)-END2(3)
          B=EYE(3)*(EYE(2)-END2(2))
          C=END1(3)-END2(3)
          D=EYE(3)*(END1(2)-END2(2))
        ENDIF
        IF (ABS(C).LT.SMALL) C=SMALL
        D=D/C
        IF (ABS(D).LT.SMALL) D=SMALL
        DXI=XDIF/YDIF
        JY1=NINT(Y1)
        JY2=NINT(Y2)
        XI=X1+DXI*(FLOAT(JY1)-Y1)
        DYY=1.0/DY
        IF (LDASH) THEN
          DXII=(NDASH1-NDASH2-1)*DXI
          DO 50 JJ=JY1,JY2,NDASH1
            JJJ=MIN(JJ+NDASH2,JY2)
            YY=YBLC-EYE(2)+DYY*FLOAT(JJ-1)
            DO 40 J=JJ,JJJ
              K=NXP*(J-1)+NINT(XI)
              Z=ZW1+(A*YY+B)/(YY+D)
              IF (Z.GT.SBBUFF(K)) THEN
                SBBUFF(K)=Z
                SBBUFF(KSTART+K)=COLOUR
              ENDIF
              XI=XI+DXI
              YY=YY+DYY
  40        CONTINUE
            XI=XI+DXII
  50      CONTINUE
        ELSE
          YY=YBLC-EYE(2)+DYY*FLOAT(JY1-1)
          DO 60 J=JY1,JY2
            K=NXP*(J-1)+NINT(XI)
            Z=ZW1+(A*YY+B)/(YY+D)
            IF (Z.GT.SBBUFF(K)) THEN
              SBBUFF(K)=Z
              SBBUFF(KSTART+K)=COLOUR
            ENDIF
            XI=XI+DXI
            YY=YY+DYY
  60      CONTINUE
        ENDIF
      ENDIF
      END SUBROUTINE SBLINE

      SUBROUTINE SBLIN1(EYE,X,Y,Z,XW,YW)
C     ----------------------------------
C
      REAL EYE(*)
C
      XLAM=EYE(3)/(EYE(3)-Z)
      XW=EYE(1)+XLAM*(X-EYE(1))
      YW=EYE(2)+XLAM*(Y-EYE(2))
      END SUBROUTINE SBLIN1

      SUBROUTINE SBPLAN(EYE,NV,VERT,IC1,IC2,LIGHT)
C     --------------------------------------------
C
      REAL             EYE(*),VERT(3,*),LIGHT(*)
C
      REAL*8           XLNORM,ZZ,DZZ
      REAL             XW(400),YW(400)
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots a diffusively-lit coloured plane; the user 
C    must ensure that all the verticies lie in a flat plane, and that 
C    the bounding polygon be convex (so that the angle at any vertex
C    <= 180 degs). All (x,y,z) values are taken to be given in world 
C    coordinates. The z-component of the eye-poisition should be 
C    positive and that of the vertices should be negative; the viewing-
C    screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    NV       R*4    I       -      No. of verticies (>=3).
C    VERT     R*4    I     3 x NV   (x,y,z) coordinate of verticies.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBLIN1     Calculates the projection of (x,y,z) on viewing screen.
C
C History
C   D. S. Sivia       4 Apr 1995  Initial release.
C   D. S. Sivia      24 Oct 1997 "Safe-guarded" some rounding errors.
C-----------------------------------------------------------------------
C
C Carry out some initial checks and calculate the coordinates of the 
C projected polygon.
C
      SMALL=1.0E-10
      SMALL2=SMALL**2
      IF (EYE(3).LE.SMALL) RETURN
      IF (NV.LT.3 .OR. NV.GT.400) RETURN
      DO 10 I=1,NV
  10    IF (VERT(3,I).GE.0.0) RETURN
      XMIN=+1.0E20
      XMAX=-1.0E20
      YMIN=+1.0E20
      YMAX=-1.0E20
      DO 20 I=1,NV
        CALL SBLIN1(EYE,VERT(1,I),VERT(2,I),VERT(3,I),XW(I),YW(I))
        IF (XW(I).LT.XMIN) THEN
          XMIN=XW(I)
          ILEFT=I
        ENDIF
        IF (YW(I).LT.YMIN) THEN
          YMIN=YW(I)
          JBOTOM=I
        ENDIF
        XMAX=MAX(XW(I),XMAX)
        YMAX=MAX(YW(I),YMAX)
  20  CONTINUE
      IF (XMIN.GE.XTRC .OR. XMAX.LE.XBLC) RETURN
      IF (YMIN.GE.YTRC .OR. YMAX.LE.YBLC) RETURN
C
C Find the outward normal seen by the eye, and activate the appropriate 
C colour.
C
      AX=VERT(1,2)-VERT(1,1)
      AY=VERT(2,2)-VERT(2,1)
      AZ=VERT(3,2)-VERT(3,1)
      BX=VERT(1,1)-VERT(1,NV)
      BY=VERT(2,1)-VERT(2,NV)
      BZ=VERT(3,1)-VERT(3,NV)
      XN=BY*AZ-AY*BZ
      YN=BZ*AX-AZ*BX
      ZN=BX*AY-AX*BY
      TEN=XN*(EYE(1)-VERT(1,1))+YN*(EYE(2)-VERT(2,1))
     *   +ZN*(EYE(3)-VERT(3,1))
      COLOUR=FLOAT(IC1)
      NC=IC2-IC1
      IF (NC.GT.0) THEN
        TNL=XN*LIGHT(1)+YN*LIGHT(2)+ZN*LIGHT(3)
        IF (TEN.LT.0.0) TNL=-TNL
        COSDIF=0.0
        IF (TNL.LT.0.0) THEN
          TN2=XN**2+YN**2+ZN**2
          TL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
          COSDIF=MIN(-TNL/SQRT(TN2*TL2+SMALL2),1.0)
        ENDIF
        COLOUR=COLOUR+COSDIF*FLOAT(NC)
      ENDIF
C
C Plot the projected polygon.
C
      XLNORM=DBLE(EYE(3))*DBLE(TEN)
      EYENRM=XN*EYE(1)+YN*EYE(2)+ZN*EYE(3)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      DYJ=1.0/DY
      DXI=1.0/DX
      SAFER=0.0001
      IF ((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
        JYMIN=INT((YMIN-YBLC)*DY)+2
        JYMAX=MIN(INT((YMAX-YBLC)*DY)+1,NYP)
        IF (JYMIN.GT.JYMAX) RETURN
        YJ=YBLC+(FLOAT(JYMIN-1)+SAFER)*DYJ
        NVL2=JBOTOM
        NVR2=JBOTOM
        J1=JYMIN
        DO 50 IVERT=1,NV
          IF (YJ.GT.YW(NVL2)) THEN
   1        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVL2)) GOTO 1
            YDIFL=YW(NVL2)-YW(NVL1)
            IF (ABS(YDIFL).LT.SMALL) YDIFL=SMALL
            GRADL=(XW(NVL2)-XW(NVL1))/YDIFL
          ENDIF
          IF (YJ.GT.YW(NVR2)) THEN
   2        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVR2)) GOTO 2
            YDIFR=YW(NVR2)-YW(NVR1)
            IF (ABS(YDIFR).LT.SMALL) YDIFR=SMALL
            GRADR=(XW(NVR2)-XW(NVR1))/YDIFR
          ENDIF
          IF (YW(NVL2).LT.YW(NVR2)) THEN
            J2=MIN(INT((YW(NVL2)-YBLC)*DY)+1,JYMAX)
          ELSE
            J2=MIN(INT((YW(NVR2)-YBLC)*DY)+1,JYMAX)
          ENDIF
          DO 40 J=J1,J2
            IF (J.GE.1) THEN
              XL=XW(NVL1)+GRADL*(YJ-YW(NVL1))
              XR=XW(NVR1)+GRADR*(YJ-YW(NVR1))
              ISTEP=1
              IX1=MAX(INT((XL-XBLC)*DX)+2,1)
              IX2=MIN(INT((XR-XBLC)*DX)+1,NXP)
              IF (IX1.GT.IX2) THEN
                ISTEP=-1
                IX1=MIN(IX1-1,NXP)
                IX2=MAX(IX2+1,1)
              ENDIF
              DZZ=DBLE(FLOAT(ISTEP)*DXI*XN)
              ZZ=DBLE(EYENRM-(XBLC+FLOAT(IX1-1)*DXI)*XN-YJ*YN)
              K=(J-1)*NXP+IX1
              DO 30 I=IX1,IX2,ISTEP
                Z=EYE(3)-SNGL(XLNORM/ZZ)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  SBBUFF(KSTART+K)=COLOUR
                ENDIF
                ZZ=ZZ-DZZ
                K=K+ISTEP
  30          CONTINUE
            ENDIF
            YJ=YJ+DYJ
  40      CONTINUE
          J1=J2+1
          IF (J1.GT.JYMAX) RETURN
  50    CONTINUE
      ELSE
        IXMIN=INT((XMIN-XBLC)*DX)+2
        IXMAX=MIN(INT((XMAX-XBLC)*DX)+1,NXP)
        IF (IXMIN.GT.IXMAX) RETURN
        XI=XBLC+(FLOAT(IXMIN-1)+SAFER)*DXI
        NVL2=ILEFT
        NVR2=ILEFT
        I1=IXMIN
        DO 80 IVERT=1,NV
          IF (XI.GT.XW(NVL2)) THEN
   3        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVL2)) GOTO 3
            XDIFL=XW(NVL2)-XW(NVL1)
            IF (ABS(XDIFL).LT.SMALL) XDIFL=SMALL
            GRADL=(YW(NVL2)-YW(NVL1))/XDIFL
          ENDIF
          IF (XI.GT.XW(NVR2)) THEN
   4        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVR2)) GOTO 4
            XDIFR=XW(NVR2)-XW(NVR1)
            IF (ABS(XDIFR).LT.SMALL) XDIFR=SMALL
            GRADR=(YW(NVR2)-YW(NVR1))/XDIFR
          ENDIF
          IF (XW(NVL2).LT.XW(NVR2)) THEN
            I2=MIN(INT((XW(NVL2)-XBLC)*DX)+1,IXMAX)
          ELSE
            I2=MIN(INT((XW(NVR2)-XBLC)*DX)+1,IXMAX)
          ENDIF
          DO 70 I=I1,I2
            IF (I.GE.1) THEN
              YL=YW(NVL1)+GRADL*(XI-XW(NVL1))
              YR=YW(NVR1)+GRADR*(XI-XW(NVR1))
              ISTEP=1
              JY1=MAX(INT((YL-YBLC)*DY)+2,1)
              JY2=MIN(INT((YR-YBLC)*DY)+1,NYP)
              IF (JY1.GT.JY2) THEN
                ISTEP=-1
                JY1=MIN(JY1-1,NYP)
                JY2=MAX(JY2+1,1)
              ENDIF
              DZZ=DBLE(FLOAT(ISTEP)*DYJ*YN)
              ZZ=DBLE(EYENRM-(YBLC+FLOAT(JY1-1)*DYJ)*YN-XI*XN)
              K=(JY1-1)*NXP+I
              KSTEP=ISTEP*NXP
              DO 60 J=JY1,JY2,ISTEP
                Z=EYE(3)-SNGL(XLNORM/ZZ)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  SBBUFF(KSTART+K)=COLOUR
                ENDIF
                ZZ=ZZ-DZZ
                K=K+KSTEP
  60          CONTINUE
            ENDIF
            XI=XI+DXI
  70      CONTINUE
          I1=I2+1
          IF (I1.GT.IXMAX) RETURN
  80    CONTINUE
      ENDIF
      END SUBROUTINE SBPLAN

      SUBROUTINE SBROD(EYE,END1,END2,RADIUS,IC1,IC2,LIGHT,NSIDES,LEND)
C     ----------------------------------------------------------------
C
      REAL     EYE(*),END1(*),END2(*),LIGHT(*)
      LOGICAL  LEND
C
      LOGICAL  LEND1,LSTART
      REAL     ENDVRT(3,361),SIDVRT(3,4)
      REAL     SINROT(361),COSROT(361)
      SAVE     LSTART,NROT0,SINROT,COSROT,SCLNRM
      DATA     NRTMAX,PI,LSTART /361,3.141592654,.FALSE./
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
!




C      This subroutine plots a diffusively-shaded coloured rod. All 
C    (x,y,z) values are taken to be given in world coordinates. The 
C    z-component of the eye-poisition should be positive and that of
C    the rod-ends should be negative (< -radius); the viewing-screen
C    is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    END1     R*4    I       3      (x,y,z) coordinate of rod-end 1.
C    END2     R*4    I       3      (x,y,z) coordinate of rod-end 2.
C    RADIUS   R*4    I       -      Radius of cylinderical rod.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    NSIDES   I*4    I       -      The order of the polygon to be used
C                                   for the cross-section of the rod.
C    LEND     L*1    I       -      If true, plot the end of the rod.
C
C External Calls
C     SBLINE     Draws a straight line between two points.
C     SBPLAN     Plots a coloured plane.
C     SBROD1     Initialises the array of sines and coses, if necessary.
C     SBRCOP     Copies one real array to another.
C
C History
C   D. S. Sivia       4 Apr 1995  Initial release.
C   D. S. Sivia      22 Oct 1997  Increased NRTMAX from 73 to 361.
C-----------------------------------------------------------------------
C
C Some initial checks.
C
      IF (NSIDES.LE.2) THEN
        CALL SBLINE(EYE,END1,END2,IC2,.FALSE.)
        RETURN
      ENDIF
      SMALL=1.0E-20
      IF (EYE(3).LE.0.0) RETURN
      IF (RADIUS.LE.0.0) RETURN
      IF (END1(3).GT.-RADIUS .OR. END2(3).GT.-RADIUS) RETURN
      XL=END2(1)-END1(1)
      YL=END2(2)-END1(2)
      ZL=END2(3)-END1(3)
      RLENG2=XL**2+YL**2+ZL**2
      IF (RLENG2.LT.SMALL) RETURN
C
C Find out which end of the rod, if any, can be seen.
C
      LEND1=.FALSE.
      X0=END1(1)
      Y0=END1(2)
      Z0=END1(3)
      EYEND=XL*(EYE(1)-X0)+YL*(EYE(2)-Y0)+ZL*(EYE(3)-Z0)
      IF (EYEND.LT.0.0) THEN
        LEND1=.TRUE.
      ELSE
        X0=END2(1)
        Y0=END2(2)
        Z0=END2(3)
        XL=-XL
        YL=-YL
        ZL=-ZL
        EYEND=XL*(EYE(1)-X0)+YL*(EYE(2)-Y0)+ZL*(EYE(3)-Z0)
        IF (EYEND.LT.0.0) LEND1=.TRUE.
      ENDIF
      SINTHT=0.0
      COSTHT=1.0
      SINPHI=0.0
      COSPHI=1.0
      RXY2=XL**2+YL**2
      IF (RXY2.GT.SMALL) THEN
        RLENG=SQRT(RLENG2)
        RXY=SQRT(RXY2)
        SINTHT=+RXY/RLENG
        COSTHT=-ZL/RLENG  
        SINPHI=-YL/RXY
        COSPHI=-XL/RXY
      ENDIF
      IF (.NOT. LEND) LEND1=.FALSE.
C
C Sweep around the rod and plot the shaded surface.
C
      NROT=MIN(NSIDES+1,NRTMAX)
      CALL SBROD1(LSTART,NROT0,NROT,SINROT,COSROT,PI,SCLNRM)
      DO 10 I=1,NROT
        X=RADIUS*COSROT(I)
        Y=RADIUS*SINROT(I)
        ENDVRT(1,I)=X0+X*COSTHT*COSPHI-Y*SINPHI
        ENDVRT(2,I)=Y0+X*COSTHT*SINPHI+Y*COSPHI
        ENDVRT(3,I)=Z0-X*SINTHT
  10  CONTINUE
      DO 20 J=2,NROT
        I=J-1
        XN=SCLNRM*(ENDVRT(1,I)+ENDVRT(1,J)-X0-X0)
        YN=SCLNRM*(ENDVRT(2,I)+ENDVRT(2,J)-Y0-Y0)
        ZN=SCLNRM*(ENDVRT(3,I)+ENDVRT(3,J)-Z0-Z0)
        ENN2=(EYE(1)-X0-XN)*XN+(EYE(2)-Y0-YN)*YN+(EYE(3)-Z0-ZN)*ZN
        IF (ENN2.GT.0.0) THEN
          CALL SBRCOP(ENDVRT(1,I),SIDVRT(1,1),3)
          SIDVRT(1,2)=ENDVRT(1,I)+XL
          SIDVRT(2,2)=ENDVRT(2,I)+YL
          SIDVRT(3,2)=ENDVRT(3,I)+ZL
          SIDVRT(1,3)=ENDVRT(1,J)+XL
          SIDVRT(2,3)=ENDVRT(2,J)+YL
          SIDVRT(3,3)=ENDVRT(3,J)+ZL
          CALL SBRCOP(ENDVRT(1,J),SIDVRT(1,4),3)
          CALL SBPLAN(EYE,4,SIDVRT,IC1,IC2,LIGHT)
        ENDIF
  20  CONTINUE
      IF (LEND1) CALL SBPLAN(EYE,NROT-1,ENDVRT,IC1,IC2,LIGHT)
      END SUBROUTINE SBROD

      SUBROUTINE SBROD1(LSTART,NROT0,NROT,SINROT,COSROT,PI,SCLNRM)
C     ------------------------------------------------------------
C
      REAL    COSROT(*),SINROT(*)
      LOGICAL LSTART
C
      IF (LSTART .AND. NROT.EQ.NROT0) RETURN
      ROT=0.0
      DROT=2.0*PI/FLOAT(NROT-1)
      DO 10 I=1,NROT
        SINROT(I)=SIN(ROT)
        COSROT(I)=COS(ROT)
        ROT=ROT+DROT
  10  CONTINUE
      SCLNRM=0.5/COS(DROT/2.0)
      LSTART=.TRUE.
      NROT0=NROT
      END SUBROUTINE SBROD1

      SUBROUTINE SBCONE(EYE,BASE,APEX,RADIUS,IC1,IC2,LIGHT,NSIDES)
C     ------------------------------------------------------------
C
      REAL     EYE(*),BASE(*),APEX(*),LIGHT(*)
C
      LOGICAL  LBASE,LSTART
      REAL     BASVRT(3,361),SIDVRT(3,3)
      REAL     SINROT(361),COSROT(361)
      SAVE     LSTART,NROT0,SINROT,COSROT,SCLNRM
      DATA     NRTMAX,PI,LSTART /361,3.141592654,.FALSE./
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots a diffusively-shaded coloured right-angular
C    cone. All (x,y,z) values are taken to be given in world coordinates.
C    The z-component of the eye-poisition should be positive and that of
C    the base and appex of the cone should be negative (< -radius); the 
C    viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    BASE     R*4    I       3      (x,y,z) coordinate of the centre of
C                                   the base of the cone.
C    APEX     R*4    I       3      (x,y,z) coordinate of the apex.
C    RADIUS   R*4    I       -      Radius of the base of the cone.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    NSIDES   I*4    I       -      The order of the polygon to be used
C                                   for the cross-section of the cone.
C
C External Calls
C     SBLINE     Draws a straight line between two points.
C     SBPLAN     Plots a coloured plane.
C     SBROD1     Initialises the array of sines and coses, if necessary.
C     SBRCOP     Copies one real array to another.
C
C History
C   D. S. Sivia      29 Jun 1995  Initial release.
C   D. S. Sivia      22 Oct 1997  Increased NRTMAX from 73 to 361.
C-----------------------------------------------------------------------
C
C Some initial checks.
C
      SMALL=1.0E-20
      IF (NSIDES.LE.2) RETURN
      IF (EYE(3).LE.0.0) RETURN
      IF (RADIUS.LE.0.0) RETURN
      IF (BASE(3).GT.-RADIUS .OR. APEX(3).GE.0.0) RETURN
      XL=APEX(1)-BASE(1)
      YL=APEX(2)-BASE(2)
      ZL=APEX(3)-BASE(3)
      RLENG2=XL**2+YL**2+ZL**2
      IF (RLENG2.LT.SMALL) RETURN
C
C Find out whether the base of the cone can be seen.
C
      LBASE=.FALSE.
      X0=BASE(1)
      Y0=BASE(2)
      Z0=BASE(3)
      EYEND=XL*(EYE(1)-X0)+YL*(EYE(2)-Y0)+ZL*(EYE(3)-Z0)
      IF (EYEND.LT.0.0) LBASE=.TRUE.
      SINTHT=0.0
      COSTHT=1.0
      SINPHI=0.0
      COSPHI=1.0
      RXY2=XL**2+YL**2
      IF (RXY2.GT.SMALL) THEN
        RLENG=SQRT(RLENG2)
        RXY=SQRT(RXY2)
        SINTHT=+RXY/RLENG
        COSTHT=-ZL/RLENG  
        SINPHI=-YL/RXY
        COSPHI=-XL/RXY
      ENDIF
C
C Sweep around the rod and plot the shaded surface.
C
      NROT=MIN(NSIDES+1,NRTMAX)
      CALL SBROD1(LSTART,NROT0,NROT,SINROT,COSROT,PI,SCLNRM)
      DO 10 I=1,NROT
        X=RADIUS*COSROT(I)
        Y=RADIUS*SINROT(I)
        BASVRT(1,I)=X0+X*COSTHT*COSPHI-Y*SINPHI
        BASVRT(2,I)=Y0+X*COSTHT*SINPHI+Y*COSPHI
        BASVRT(3,I)=Z0-X*SINTHT
  10  CONTINUE
      IF (LBASE) CALL SBPLAN(EYE,NROT-1,BASVRT,IC1,IC2,LIGHT)
      CALL SBRCOP(APEX,SIDVRT(1,3),3)
      DO 20 J=1,NROT-1
        CALL SBRCOP(BASVRT(1,J),SIDVRT(1,1),6)
        CALL SBPLAN(EYE,3,SIDVRT,IC1,IC2,LIGHT)
  20  CONTINUE
      END SUBROUTINE SBCONE

      SUBROUTINE SBSLIC(EYE,LATICE,DENS,N1,N2,N3,DLOW,DHIGH,IC1,IC2,
     *                  SLNORM,APOINT,ICEDGE)
C     --------------------------------------------------------------
C
      REAL             EYE(*),LATICE(3,*),DENS(0:N1,0:N2,0:N3)
      REAL             SLNORM(*),APOINT(*)
C
      REAL*8           XLNORM,ZZ,DZZ
      REAL             BAS(3,3),MTRX(3,3)
      REAL             END1(3),END2(3),VERT(3,12),XW(20),YW(20)
      LOGICAL          LPS,LCOLOR,LVERT(12)
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
!




C      This subroutine plots a "grey-scale" slice through a unit-cell
C    of density. All (x,y,z) values are taken to be given in world 
C    coordinates. The z-component of the eye-poisition should be 
C    positive and that of all the lattice-vertices should be negative; 
C    the viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    LATICE   R*4    I     3 x 4    (x,y,z) coordinates of the origin
C                                   and the a, b & C lattice-vertices.
C    DENS     R*4    I     (N1+1)   The density at regular points within
C                        x (N2+1)   the unit cell, wrapped around so
C                        x (N3+1)   that DENS(0,J,K)=DENS(N1,J,K) etc..
C    N1,N2,N3 I*4    I       -      The dimensions of the unit-cell grid.
C    DLOW     R*4    I       -      Density for the lowest colour-index.
C    DHIGH    R*4    I       -      Density for the highest colour-index.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    SLNORM   R*4    I       3      (x,y,z) direction of the normal to 
C                                   the slice to be "grey-scaled".
C    APONIT   R*4    I       3      (x,y,z) coordinate of a point within
C                                   the slice to be "grey-scaled".
C    ICEDGE   I*4    I       -      If >=0, it's the colour-index for the
C                                   boundary of the "grey-scaled" slice.
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBSLC1     Plots a side of the unit cell and sees if it is cut by
C                the slice to be "grey-scaled".
C     SBSLC2     Calculates the coordinates of the projected polygon.
C     SLSLC3     Calculates the appropriate colour for a given pixel.
C     SBLIN1     Calculates the projection of (x,y,z) on viewing screen.
C
C History
C   D. S. Sivia      30 Apr 1995  Initial release.
C   D. S. Sivia       7 Jul 1995  Fixed some bug.
C   D. S. Sivia      24 Oct 1997 "Safe-guarded" some rounding errors.
C-----------------------------------------------------------------------
C
C Carry out some initial checks.
C
      SMALL=1.0E-10
      SMALL2=SMALL**2
      IF (EYE(3).LE.SMALL) RETURN
      IF (N1.LT.1 .OR. N2.LT.1 .OR. N3.LT.1) RETURN
      SNRM=1.0/SQRT(SLNORM(1)**2+SLNORM(2)**2+SLNORM(3)**2+SMALL2)
      XN=SLNORM(1)*SNRM
      YN=SLNORM(2)*SNRM
      ZN=SLNORM(3)*SNRM
      XAE=EYE(1)-APOINT(1)
      YAE=EYE(2)-APOINT(2)
      ZAE=EYE(3)-APOINT(3)
      XLNORM=DBLE(XN*XAE)+DBLE(YN*YAE)+DBLE(ZN*ZAE)
      COSNRM=SNGL(XLNORM/SQRT(DBLE(XAE**2)+DBLE(YAE**2)+DBLE(ZAE**2)
     *      +DBLE(SMALL2)))
      IF (ABS(COSNRM).LT.0.001) RETURN
      DO 10 J=1,3
        BAS(1,J)=LATICE(1,J+1)-LATICE(1,1)
        BAS(2,J)=LATICE(2,J+1)-LATICE(2,1)
        BAS(3,J)=LATICE(3,J+1)-LATICE(3,1)
        BAS2J=BAS(1,J)**2+BAS(2,J)**2+BAS(3,J)**2
        IF (BAS2J.LT.SMALL2) RETURN
  10  CONTINUE
      COL0=FLOAT(IC1)
      CSCL=FLOAT(IC2-IC1)
      COLNRM=DHIGH-DLOW
      IF (ABS(COLNRM).LT.SMALL) RETURN
      DCOL=1.0/COLNRM
C
C Set up matrix for real-space to lattice-index transformation.
C
      DET=BAS(1,1)*BAS(2,2)*BAS(3,3)+BAS(1,2)*BAS(2,3)*BAS(3,1)
     *   +BAS(1,3)*BAS(2,1)*BAS(3,2)-BAS(3,1)*BAS(2,2)*BAS(1,3)
     *   -BAS(3,2)*BAS(2,3)*BAS(1,1)-BAS(3,3)*BAS(2,1)*BAS(1,2)
      IF (ABS(DET).LT.SMALL2) RETURN
      DETNRM=1.0/DET
      MTRX(1,1)=DETNRM*(BAS(2,2)*BAS(3,3)-BAS(2,3)*BAS(3,2))
      MTRX(1,2)=DETNRM*(BAS(2,3)*BAS(3,1)-BAS(2,1)*BAS(3,3))
      MTRX(1,3)=DETNRM*(BAS(2,1)*BAS(3,2)-BAS(2,2)*BAS(3,1))
      MTRX(2,1)=DETNRM*(BAS(3,2)*BAS(1,3)-BAS(3,3)*BAS(1,2))
      MTRX(2,2)=DETNRM*(BAS(3,3)*BAS(1,1)-BAS(3,1)*BAS(1,3))
      MTRX(2,3)=DETNRM*(BAS(3,1)*BAS(1,2)-BAS(3,2)*BAS(1,1))
      MTRX(3,1)=DETNRM*(BAS(1,2)*BAS(2,3)-BAS(1,3)*BAS(2,2))
      MTRX(3,2)=DETNRM*(BAS(1,3)*BAS(2,1)-BAS(1,1)*BAS(2,3))
      MTRX(3,3)=DETNRM*(BAS(1,1)*BAS(2,2)-BAS(1,2)*BAS(2,1))
C
C Draw the frame of the unit cell and calculate the coordinates of the 
C projected polygon.
C
      NVERT=0
      II=0
      DO 20 L=1,12,4
        I=II+2
        J=MOD(II+1,3)+2
        K=MOD(II+2,3)+2
        CALL SBSLC1(LATICE,LATICE(1,K),XN,YN,ZN,APOINT,LVERT(L),
     *              VERT(1,L),NVERT)
        END1(1)=LATICE(1,I)+LATICE(1,J)-LATICE(1,1)
        END1(2)=LATICE(2,I)+LATICE(2,J)-LATICE(2,1)
        END1(3)=LATICE(3,I)+LATICE(3,J)-LATICE(3,1)
        CALL SBSLC1(LATICE(1,I),END1,XN,YN,ZN,APOINT,LVERT(L+1),
     *              VERT(1,L+1),NVERT)
        CALL SBSLC1(LATICE(1,J),END1,XN,YN,ZN,APOINT,LVERT(L+2),
     *              VERT(1,L+2),NVERT)
        END2(1)=END1(1)+LATICE(1,K)-LATICE(1,1)
        END2(2)=END1(2)+LATICE(2,K)-LATICE(2,1)
        END2(3)=END1(3)+LATICE(3,K)-LATICE(3,1)
        CALL SBSLC1(END1,END2,XN,YN,ZN,APOINT,LVERT(L+3),
     *              VERT(1,L+3),NVERT)
        II=II+1
  20  CONTINUE
      IF (NVERT.LT.3) RETURN
      CALL SBSLC2(EYE,LVERT,VERT,NVERT,XN,YN,ZN,XW,YW,ICEDGE,ZDLINE)
C
C Paint the projected polygon slice.
C
      XMIN=+1.0E20
      XMAX=-1.0E20
      YMIN=+1.0E20
      YMAX=-1.0E20
      DO 30 I=1,NVERT
        IF (XW(I).LT.XMIN) THEN
          XMIN=XW(I)
          ILEFT=I
        ENDIF
        IF (YW(I).LT.YMIN) THEN
          YMIN=YW(I)
          JBOTOM=I
        ENDIF
        XMAX=MAX(XW(I),XMAX)
        YMAX=MAX(YW(I),YMAX)
  30  CONTINUE
      IF (XMIN.GE.XTRC .OR. XMAX.LE.XBLC) RETURN
      IF (YMIN.GE.YTRC .OR. YMAX.LE.YBLC) RETURN
C
      EYENRM=XN*EYE(1)+YN*EYE(2)+ZN*EYE(3)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      DYJ=1.0/DY
      DXI=1.0/DX
      SAFER=0.0001
      IF ((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
        JYMIN=INT((YMIN-YBLC)*DY)+2
        JYMAX=MIN(INT((YMAX-YBLC)*DY)+1,NYP)
        IF (JYMIN.GT.JYMAX) RETURN
        YJ=YBLC+(FLOAT(JYMIN-1)+SAFER)*DYJ
        NVL2=JBOTOM
        NVR2=JBOTOM
        J1=JYMIN
        DO 60 IVERT=1,NVERT
          IF (YJ.GT.YW(NVL2)) THEN
   1        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NVERT
            IF (NVL2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVL2)) GOTO 1
            YDIFL=YW(NVL2)-YW(NVL1)
            IF (ABS(YDIFL).LT.SMALL) YDIFL=SMALL
            GRADL=(XW(NVL2)-XW(NVL1))/YDIFL
          ENDIF
          IF (YJ.GT.YW(NVR2)) THEN
   2        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NVERT) NVR2=1
            IF (NVR2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVR2)) GOTO 2
            YDIFR=YW(NVR2)-YW(NVR1)
            IF (ABS(YDIFR).LT.SMALL) YDIFR=SMALL
            GRADR=(XW(NVR2)-XW(NVR1))/YDIFR
          ENDIF
          IF (YW(NVL2).LT.YW(NVR2)) THEN
            J2=MIN(INT((YW(NVL2)-YBLC)*DY)+1,JYMAX)
          ELSE
            J2=MIN(INT((YW(NVR2)-YBLC)*DY)+1,JYMAX)
          ENDIF
          DO 50 J=J1,J2
            IF (J.GE.1) THEN
              XL=XW(NVL1)+GRADL*(YJ-YW(NVL1))
              XR=XW(NVR1)+GRADR*(YJ-YW(NVR1))
              ISTEP=1
              IX1=MAX(INT((XL-XBLC)*DX)+2,1)
              IX2=MIN(INT((XR-XBLC)*DX)+1,NXP)
              IF (IX1.GT.IX2) THEN
                ISTEP=-1
                IX1=MIN(IX1-1,NXP)
                IX2=MAX(IX2+1,1)
              ENDIF
              XI=XBLC+FLOAT(IX1-1)*DXI
              SDXI=FLOAT(ISTEP)*DXI
              DZZ=DBLE(SDXI*XN)
              ZZ=DBLE(EYENRM-XI*XN-YJ*YN)
              K=(J-1)*NXP+IX1
              DO 40 I=IX1,IX2,ISTEP
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF ((Z-SBBUFF(K)).GT.ZDLINE) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))-LATICE(1,1)
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))-LATICE(2,1)
                  Z=Z-LATICE(3,1)
                  CALL SBSLC3(DENS,N1,N2,N3,X,Y,Z,MTRX,DLOW,DCOL,COLOUR)
                  SBBUFF(KSTART+K)=COL0+CSCL*COLOUR
                ENDIF
                XI=XI+SDXI
                ZZ=ZZ-DZZ
                K=K+ISTEP
  40          CONTINUE
            ENDIF
            YJ=YJ+DYJ
  50      CONTINUE
          J1=J2+1
          IF (J1.GT.JYMAX) RETURN
  60    CONTINUE
      ELSE
        IXMIN=INT((XMIN-XBLC)*DX)+2
        IXMAX=MIN(INT((XMAX-XBLC)*DX)+1,NXP)
        IF (IXMIN.GT.IXMAX) RETURN
        XI=XBLC+(FLOAT(IXMIN-1)+SAFER)*DXI
        NVL2=ILEFT
        NVR2=ILEFT
        I1=IXMIN
        DO 90 IVERT=1,NVERT
          IF (XI.GT.XW(NVL2)) THEN
   3        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NVERT
            IF (NVL2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVL2)) GOTO 3
            XDIFL=XW(NVL2)-XW(NVL1)
            IF (ABS(XDIFL).LT.SMALL) XDIFL=SMALL
            GRADL=(YW(NVL2)-YW(NVL1))/XDIFL
          ENDIF
          IF (XI.GT.XW(NVR2)) THEN
   4        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NVERT) NVR2=1
            IF (NVR2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVR2)) GOTO 4
            XDIFR=XW(NVR2)-XW(NVR1)
            IF (ABS(XDIFR).LT.SMALL) XDIFR=SMALL
            GRADR=(YW(NVR2)-YW(NVR1))/XDIFR
          ENDIF
          IF (XW(NVL2).LT.XW(NVR2)) THEN
            I2=MIN(INT((XW(NVL2)-XBLC)*DX)+1,IXMAX)
          ELSE
            I2=MIN(INT((XW(NVR2)-XBLC)*DX)+1,IXMAX)
          ENDIF
          DO 80 I=I1,I2
            IF (I.GE.1) THEN
              YL=YW(NVL1)+GRADL*(XI-XW(NVL1))
              YR=YW(NVR1)+GRADR*(XI-XW(NVR1))
              ISTEP=1
              JY1=MAX(INT((YL-YBLC)*DY)+2,1)
              JY2=MIN(INT((YR-YBLC)*DY)+1,NYP)
              IF (JY1.GT.JY2) THEN
                ISTEP=-1
                JY1=MIN(JY1-1,NYP)
                JY2=MAX(JY2+1,1)
              ENDIF
              YJ=YBLC+FLOAT(JY1-1)*DYJ
              SDYJ=FLOAT(ISTEP)*DYJ
              DZZ=DBLE(SDYJ*YN)
              ZZ=DBLE(EYENRM-YJ*YN-XI*XN)
              K=(JY1-1)*NXP+I
              KSTEP=ISTEP*NXP
              DO 70 J=JY1,JY2,ISTEP
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF ((Z-SBBUFF(K)).GT.ZDLINE) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))-LATICE(1,1)
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))-LATICE(2,1)
                  Z=Z-LATICE(3,1)
                  CALL SBSLC3(DENS,N1,N2,N3,X,Y,Z,MTRX,DLOW,DCOL,COLOUR)
                  SBBUFF(KSTART+K)=COL0+CSCL*COLOUR
                ENDIF
                YJ=YJ+SDYJ
                ZZ=ZZ-DZZ
                K=K+KSTEP
  70          CONTINUE
            ENDIF
            XI=XI+DXI
  80      CONTINUE
          I1=I2+1
          IF (I1.GT.IXMAX) RETURN
  90    CONTINUE
      ENDIF
      END SUBROUTINE SBSLIC

      SUBROUTINE SBSLC1(END1,END2,XN,YN,ZN,APOINT,LVERT,VERT,NVERT)
C     -------------------------------------------------------------
C
      REAL    END1(*),END2(*),APOINT(*),VERT(*)
      LOGICAL LVERT
C
      LVERT=.FALSE.
      X12=END2(1)-END1(1)
      Y12=END2(2)-END1(2)
      Z12=END2(3)-END1(3)
      DENOM=XN*X12+YN*Y12+ZN*Z12
      COSNRM=DENOM/SQRT(X12**2+Y12**2+Z12**2+1.0E-20)
      IF (ABS(COSNRM).LT.0.001) RETURN
      XLAM=(XN*(APOINT(1)-END1(1))+YN*(APOINT(2)-END1(2))
     *     +ZN*(APOINT(3)-END1(3)))/DENOM
      IF (XLAM.GE.0.0 .AND. XLAM.LE.1.0) THEN
        LVERT=.TRUE.
        NVERT=NVERT+1
        VERT(1)=END1(1)+XLAM*X12
        VERT(2)=END1(2)+XLAM*Y12
        VERT(3)=END1(3)+XLAM*Z12
        IF (VERT(3).GE.0.0) NVERT=-1000
      ENDIF
      END SUBROUTINE SBSLC1

      SUBROUTINE SBSLC2(EYE,LVERT,VERT,NVERT,XN,YN,ZN,XW,YW,ICOL,ZDIF)
C     ----------------------------------------------------------------
C
      REAL    EYE(*),VERT(3,*),XW(*),YW(*),ANGLE(12)
      INTEGER ISORT(12)
      LOGICAL LVERT(*)
C
      IV1=0
      XBAR=0.0
      YBAR=0.0
      ZBAR=0.0
      ZMIN=+1.0E20
      ZMAX=-1.0E20
      DO 10 K=1,12
        ZMIN=MIN(ZMIN,VERT(3,K))
        ZMAX=MAX(ZMAX,VERT(3,K))
        IF (LVERT(K)) THEN
          IF (IV1.LE.0) IV1=K
          XBAR=XBAR+VERT(1,K)
          YBAR=YBAR+VERT(2,K)
          ZBAR=ZBAR+VERT(3,K)
        ENDIF
  10  CONTINUE
      ZDIF=(ZMAX-ZMIN)/5000.0
      XBAR=XBAR/FLOAT(NVERT)
      YBAR=YBAR/FLOAT(NVERT)
      ZBAR=ZBAR/FLOAT(NVERT)
      XREF=VERT(1,IV1)-XBAR
      YREF=VERT(2,IV1)-YBAR
      ZREF=VERT(3,IV1)-ZBAR
      REFNRM=1.0/SQRT(XREF**2+YREF**2+ZREF**2+1.0E-20)
      XREF=XREF*REFNRM
      YREF=YREF*REFNRM
      ZREF=ZREF*REFNRM
      XNRM=YREF*ZN-YN*ZREF
      YNRM=ZREF*XN-ZN*XREF
      ZNRM=XREF*YN-XN*YREF
      J=1
      ANGLE(J)=0.0
      ISORT(J)=IV1
      CALL SBLIN1(EYE,VERT(1,IV1),VERT(2,IV1),VERT(3,IV1),XW(J),YW(J))
      DO 40 K=IV1+1,12
        IF (LVERT(K)) THEN
          J=J+1
          XVEC=VERT(1,K)-XBAR
          YVEC=VERT(2,K)-YBAR
          ZVEC=VERT(3,K)-ZBAR
          X=XVEC*XREF+YVEC*YREF+ZVEC*ZREF
          Y=XVEC*XNRM+YVEC*YNRM+ZVEC*ZNRM
          ANGJ=ATAN2(Y,X)
          CALL SBLIN1(EYE,VERT(1,K),VERT(2,K),VERT(3,K),XWJ,YWJ)
          DO 20 I=1,J-1
  20        IF (ANGJ.LT.ANGLE(I)) GOTO 1
   1      II=I
           DO 30 I=J,II+1,-1
            XW(I)=XW(I-1)
            YW(I)=YW(I-1)
            ANGLE(I)=ANGLE(I-1)
            ISORT(I)=ISORT(I-1)
  30      CONTINUE
          XW(II)=XWJ
          YW(II)=YWJ
          ANGLE(II)=ANGJ
          ISORT(II)=K
        ENDIF
  40  CONTINUE
      IF (ICOL.GE.0.0) THEN
        DO 50 I=1,NVERT-1
          J=ISORT(I)
          K=ISORT(I+1)
          CALL SBLINE(EYE,VERT(1,J),VERT(1,K),ICOL,.FALSE.)
  50    CONTINUE
        CALL SBLINE(EYE,VERT(1,K),VERT(1,ISORT(1)),ICOL,.FALSE.)
      ENDIF
      END SUBROUTINE SBSLC2

      SUBROUTINE SBSLC3(DENS,N1,N2,N3,X,Y,Z,BAS,DLOW,DCOL,COLOUR)
C     -----------------------------------------------------------
C
      REAL DENS(0:N1,0:N2,0:N3),BAS(3,*)
      DATA RMIN,RMAX /0.00001,0.99999/
C
      XI=MIN(MAX(X*BAS(1,1)+Y*BAS(2,1)+Z*BAS(3,1),RMIN),RMAX)*FLOAT(N1)
      YJ=MIN(MAX(X*BAS(1,2)+Y*BAS(2,2)+Z*BAS(3,2),RMIN),RMAX)*FLOAT(N2)
      ZK=MIN(MAX(X*BAS(1,3)+Y*BAS(2,3)+Z*BAS(3,3),RMIN),RMAX)*FLOAT(N3)
      I=INT(XI)
      J=INT(YJ)
      K=INT(ZK)
      II=I+1
      JJ=J+1
      KK=K+1
      DX=XI-FLOAT(I)
      DY=YJ-FLOAT(J)
      DZ=ZK-FLOAT(K)
      D1=(1.0-DX)*(DENS(I,J,K)+DY*(DENS(I,JJ,K)-DENS(I,J,K)))
     *       +DX*(DENS(II,J,K)+DY*(DENS(II,JJ,K)-DENS(II,J,K)))
      D2=(1.0-DX)*(DENS(I,J,KK)+DY*(DENS(I,JJ,KK)-DENS(I,J,KK)))
     *       +DX*(DENS(II,J,KK)+DY*(DENS(II,JJ,KK)-DENS(II,J,KK)))
      COLOUR=MIN(MAX((D1+DZ*(D2-D1)-DLOW)*DCOL,RMIN),RMAX)
      END SUBROUTINE SBSLC3

      SUBROUTINE SBSURF(EYE,LATICE,DENS,N1,N2,N3,DSURF,IC1,IC2,LIGHT,
     *                  LSHINE)
C     ---------------------------------------------------------------
C
      REAL             EYE(*),LATICE(3,*),DENS(0:N1,0:N2,0:N3),LIGHT(*)
      LOGICAL          LSHINE
C
      REAL*8           XLNORM,PNEYE(3),DSMAL2
      REAL             BAS(3,3),MTRX
      REAL             XYZ(3),DXYZ(3,3),FRCXYZ(12),DDXYZ(3,12,2)
      REAL             DLOCAL(8),VERT(3,12),GRDSCL(3)
      LOGICAL          LPS,LCOLOR,LEMPTY
      INTEGER          IVERT(8)
      COMMON  /SRFCOM/ GRDCUB(3,8),MTRX(3,3),ORIG(3),XL2,COL0,COLSCL
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots an iso-surface through a unit-cell of
C    density. All (x,y,z) values are taken to be given in world 
C    coordinates. The z-component of the eye-poisition should be 
C    positive and that of all the lattice-vertices should be negative; 
C    the viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    LATICE   R*4    I     3 x 4    (x,y,z) coordinates of the origin
C                                   and the a, b & C lattice-vertices.
C    DENS     R*4    I     (N1+1)   The density at regular points within
C                        x (N2+1)   the unit cell, wrapped around so
C                        x (N3+1)   that DENS(0,J,K)=DENS(N1,J,K) etc..
C    N1,N2,N3 I*4    I       -      The dimensions of the unit-cell grid.
C    DSURF    R*4    I       -      Density for the iso-surface.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    LSHINE   L*1    I       -      Shiny surface if TRUE, else diffuse.
C
C Globals 
C    SFTBUF
C    SRFCOM
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBSRF0     A quick check in case there are no iso-surafces.
C     SBSRF1     Anlayses a 2-d box in the surface of the unit-cell.
C     SBSRF2     Paints a 2-d box in the surface of the unit-cell.
C     SBSRF3     Analyses a 3-d box within the unit-cell.
C     SBSRF4     Initialises the gradients for a 3-d box.
C     SBSRF5     Breaks up the iso-surface in a 3-d box into triangles.
C     SBSRF6     Paints a triangular patch of an iso-surface.
C
C History
C   D. S. Sivia       3 May 1995  Initial release.
C   D. S. Sivia       7 Jul 1995  Fixed bug in determinant calculation.
C   D. S. Sivia      20 Oct 1995  Speeded up computations slightly.
C   D. S. Sivia      13 Dec 1995  A bit more tinkering for speed.
C   D. S. Sivia      14 Jun 1996  Completely new algorithm!
C   D. S. Sivia      24 Oct 1997 "Safe-guarded" some rounding errors.
C-----------------------------------------------------------------------
C
C Carry out some initial checks.
C
      SMALL=1.0E-10
      SMALL2=SMALL**2
      DSMAL2=DBLE(SMALL2)
      IF (EYE(3).LE.SMALL) RETURN
      IF (N1.LT.1 .OR. N2.LT.1 .OR. N3.LT.1) RETURN
      XL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
      IF (XL2.LT.SMALL) RETURN
      IF (LATICE(3,1).GE.0.0) RETURN
      ZFAR=LATICE(3,2)+LATICE(3,3)+LATICE(3,4)-2.0*LATICE(3,1)
      IF (ZFAR.GE.0.0) RETURN
      DO 10 J=1,3
        IF (LATICE(3,J+1).GE.0.0) RETURN
        BAS(1,J)=LATICE(1,J+1)-LATICE(1,1)
        BAS(2,J)=LATICE(2,J+1)-LATICE(2,1)
        BAS(3,J)=LATICE(3,J+1)-LATICE(3,1)
        IF ((ZFAR-BAS(3,J)).GE.0.0) RETURN
        BAS2J=BAS(1,J)**2+BAS(2,J)**2+BAS(3,J)**2
        IF (BAS2J.LT.SMALL2) RETURN
  10  CONTINUE
      NTOT=(N1+1)*(N2+1)*(N3+1)
      CALL SBSRF0(DENS,NTOT,DSURF,LEMPTY)
      IF (LEMPTY) RETURN
C
C Set up matrix for real-space to lattice-index transformation.
C
      XN1=0.99999*FLOAT(N1)
      XN2=0.99999*FLOAT(N2)
      XN3=0.99999*FLOAT(N3)
      DET=BAS(1,1)*BAS(2,2)*BAS(3,3)+BAS(1,2)*BAS(2,3)*BAS(3,1)
     *   +BAS(1,3)*BAS(2,1)*BAS(3,2)-BAS(3,1)*BAS(2,2)*BAS(1,3)
     *   -BAS(3,2)*BAS(2,3)*BAS(1,1)-BAS(3,3)*BAS(2,1)*BAS(1,2)
      IF (ABS(DET).LT.SMALL2) RETURN
      DETNRM=1.0/DET
      MTRX(1,1)=XN1*DETNRM*(BAS(2,2)*BAS(3,3)-BAS(2,3)*BAS(3,2))
      MTRX(1,2)=XN2*DETNRM*(BAS(2,3)*BAS(3,1)-BAS(2,1)*BAS(3,3))
      MTRX(1,3)=XN3*DETNRM*(BAS(2,1)*BAS(3,2)-BAS(2,2)*BAS(3,1))
      MTRX(2,1)=XN1*DETNRM*(BAS(3,2)*BAS(1,3)-BAS(3,3)*BAS(1,2))
      MTRX(2,2)=XN2*DETNRM*(BAS(3,3)*BAS(1,1)-BAS(3,1)*BAS(1,3))
      MTRX(2,3)=XN3*DETNRM*(BAS(3,1)*BAS(1,2)-BAS(3,2)*BAS(1,1))
      MTRX(3,1)=XN1*DETNRM*(BAS(1,2)*BAS(2,3)-BAS(1,3)*BAS(2,2))
      MTRX(3,2)=XN2*DETNRM*(BAS(1,3)*BAS(2,1)-BAS(1,1)*BAS(2,3))
      MTRX(3,3)=XN3*DETNRM*(BAS(1,1)*BAS(2,2)-BAS(1,2)*BAS(2,1))
      CALL SBRCOP(LATICE,ORIG,3)
C
C Some general initialisations.
C
      DDSURF=MAX(ABS(DSURF),SMALL)
      IF (DSURF.LT.0.0) DDSURF=-DDSURF
      GRDSCL(1)=-0.5/(DDSURF*FLOAT(N1))
      GRDSCL(2)=-0.5/(DDSURF*FLOAT(N2))
      GRDSCL(3)=-0.5/(DDSURF*FLOAT(N3))
      COL0=FLOAT(IC1)
      COLSCL=FLOAT(IC2-IC1)
      DO 30 I=1,3
        DXYZ(I,1)=BAS(I,1)/FLOAT(N1)
        DXYZ(I,2)=BAS(I,2)/FLOAT(N2)
        DXYZ(I,3)=BAS(I,3)/FLOAT(N3)
        DDXYZ(I,1,1)=0.0
        DDXYZ(I,1,2)=DXYZ(I,1)
        DDXYZ(I,2,1)=DDXYZ(I,1,1)+DDXYZ(I,1,2)
        DDXYZ(I,2,2)=DXYZ(I,2)
        DDXYZ(I,3,1)=DDXYZ(I,2,1)+DDXYZ(I,2,2)
        DDXYZ(I,3,2)=-DXYZ(I,1)
        DDXYZ(I,4,1)=DDXYZ(I,3,1)+DDXYZ(I,3,2)
        DDXYZ(I,4,2)=-DXYZ(I,2)
        DO 20 J=1,4
          DDXYZ(I,J+4,1)=DDXYZ(I,J,1)
          DDXYZ(I,J+4,2)=DXYZ(I,3)
          DDXYZ(I,J+8,1)=DDXYZ(I,J,1)+DXYZ(I,3)
          DDXYZ(I,J+8,2)=DDXYZ(I,J,2)
  20    CONTINUE
  30  CONTINUE
C
C First paint the edges of the lattice.
C
      DO 60 IFACE=1,3
        I=IFACE
        J=MOD(IFACE,3)+1
        K=MOD(J,3)+1
        IF (IFACE.EQ.1) THEN
          IN=N1
          JN=N2
        ELSEIF (IFACE.EQ.2) THEN
          IN=N2
          JN=N3
        ELSE
          IN=N3
          JN=N1
        ENDIF
        KK=0
        XN=BAS(2,J)*BAS(3,I)-BAS(2,I)*BAS(3,J)
        YN=BAS(3,J)*BAS(1,I)-BAS(3,I)*BAS(1,J)
        ZN=BAS(1,J)*BAS(2,I)-BAS(1,I)*BAS(2,J)
        DNRM=SQRT(XN**2+YN**2+ZN**2+SMALL2)
        PNEYE(1)=DBLE(EYE(1)-0.5*(LATICE(1,I+1)+LATICE(1,J+1)))
        PNEYE(2)=DBLE(EYE(2)-0.5*(LATICE(2,I+1)+LATICE(2,J+1)))
        PNEYE(3)=DBLE(EYE(3)-0.5*(LATICE(3,I+1)+LATICE(3,J+1)))
        DEYE=SNGL(SQRT(PNEYE(1)**2+PNEYE(2)**2+PNEYE(3)**2+DSMAL2))
        XLNORM=DBLE(XN)*PNEYE(1)+DBLE(YN)*PNEYE(2)+DBLE(ZN)*PNEYE(3)
        COSSEE=SNGL(XLNORM)/(DEYE*DNRM)
        IF (COSSEE.LT.0.001) THEN
          KK=N3
          IF (IFACE.EQ.2) KK=N1
          IF (IFACE.EQ.3) KK=N2
          PNEYE(1)=PNEYE(1)+DBLE(BAS(1,K))
          PNEYE(2)=PNEYE(2)+DBLE(BAS(2,K))
          PNEYE(3)=PNEYE(3)+DBLE(BAS(3,K))
          DEYE=SNGL(SQRT(PNEYE(1)**2+PNEYE(2)**2+PNEYE(3)**2+DSMAL2))
          XLNORM=DBLE(XN)*PNEYE(1)+DBLE(YN)*PNEYE(2)+DBLE(ZN)*PNEYE(3)
          COSSEE=-SNGL(XLNORM)/(DEYE*DNRM)
        ENDIF
        IF (COSSEE.GT.0.001) THEN
          XYZ1=FLOAT(KK)*DXYZ(1,K)+LATICE(1,1)
          XYZ2=FLOAT(KK)*DXYZ(2,K)+LATICE(2,1)
          XYZ3=FLOAT(KK)*DXYZ(3,K)+LATICE(3,1)
          DO 50 J1=1,JN
            J0=J1-1
            DO 40 I1=1,IN
              I0=I1-1
              IF (IFACE.EQ.1) THEN
                DLOCAL(1)=DENS(I0,J0,KK)-DSURF
                DLOCAL(2)=DENS(I1,J0,KK)-DSURF
                DLOCAL(3)=DENS(I1,J1,KK)-DSURF
                DLOCAL(4)=DENS(I0,J1,KK)-DSURF
              ELSEIF (IFACE.EQ.2) THEN
                DLOCAL(1)=DENS(KK,I0,J0)-DSURF
                DLOCAL(2)=DENS(KK,I1,J0)-DSURF
                DLOCAL(3)=DENS(KK,I1,J1)-DSURF
                DLOCAL(4)=DENS(KK,I0,J1)-DSURF
              ELSE
                DLOCAL(1)=DENS(J0,KK,I0)-DSURF
                DLOCAL(2)=DENS(J0,KK,I1)-DSURF
                DLOCAL(3)=DENS(J1,KK,I1)-DSURF
                DLOCAL(4)=DENS(J1,KK,I0)-DSURF
              ENDIF
              CALL SBSRF1(DLOCAL,IBSIDE,FRCXYZ)
              IF (IBSIDE.NE.0) THEN
                XYZ(1)=XYZ1+DXYZ(1,I)*FLOAT(I0)
                XYZ(2)=XYZ2+DXYZ(2,I)*FLOAT(I0)
                XYZ(3)=XYZ3+DXYZ(3,I)*FLOAT(I0)
                CALL SBSRF2(XYZ,DXYZ(1,I),DXYZ(1,J),IBSIDE,FRCXYZ,VERT,
     *                      EYE,LIGHT,LSHINE)
              ENDIF
  40        CONTINUE
            XYZ1=XYZ1+DXYZ(1,J)
            XYZ2=XYZ2+DXYZ(2,J)
            XYZ3=XYZ3+DXYZ(3,J)
  50      CONTINUE
        ENDIF
  60  CONTINUE
C
C Step through each "cube" in the lattice, and paint any isosurfaces
C found therein.
C
      X00K=LATICE(1,1)
      Y00K=LATICE(2,1)
      Z00K=LATICE(3,1)
      DO 90 K1=1,N3
        K0=K1-1
        DO 80 J1=1,N2
          J0=J1-1
          X0JK=X00K+DXYZ(1,2)*FLOAT(J0)
          Y0JK=Y00K+DXYZ(2,2)*FLOAT(J0)
          Z0JK=Z00K+DXYZ(3,2)*FLOAT(J0)
          DO 70 I1=1,N1
            I0=I1-1
            DLOCAL(1)=DENS(I0,J0,K0)-DSURF
            DLOCAL(2)=DENS(I1,J0,K0)-DSURF
            DLOCAL(3)=DENS(I1,J1,K0)-DSURF
            DLOCAL(4)=DENS(I0,J1,K0)-DSURF
            DLOCAL(5)=DENS(I0,J0,K1)-DSURF
            DLOCAL(6)=DENS(I1,J0,K1)-DSURF
            DLOCAL(7)=DENS(I1,J1,K1)-DSURF
            DLOCAL(8)=DENS(I0,J1,K1)-DSURF
            CALL SBSRF3(DLOCAL,IVERT,FRCXYZ,ISUMV,ISUMF)
            IF (ISUMV.NE.0) THEN
              XYZ(1)=X0JK+DXYZ(1,1)*FLOAT(I0)
              XYZ(2)=Y0JK+DXYZ(2,1)*FLOAT(I0)
              XYZ(3)=Z0JK+DXYZ(3,1)*FLOAT(I0)
              CALL SBSRF4(DENS,N1,N2,N3,I0,J0,K0,GRDSCL,BAS,GRDCUB)
              CALL SBSRF5(XYZ,DDXYZ,ISUMV,ISUMF,IVERT,FRCXYZ,VERT,EYE,
     *                    LIGHT,LSHINE)
            ENDIF
  70      CONTINUE
  80    CONTINUE
        X00K=X00K+DXYZ(1,3)
        Y00K=Y00K+DXYZ(2,3)
        Z00K=Z00K+DXYZ(3,3)
  90  CONTINUE
      END SUBROUTINE SBSURF

      SUBROUTINE SBSRF0(DENS,NTOT,DSURF,LEMPTY)
C     -----------------------------------------
C
      REAL    DENS(*)
      LOGICAL LEMPTY
C
      LEMPTY=.TRUE.
      DO 10 I=1,NTOT
        IF (DENS(I).GT.DSURF) THEN
          LEMPTY=.FALSE.
          RETURN
        ENDIF
  10  CONTINUE
      END SUBROUTINE SBSRF0
      SUBROUTINE SBSRF1(D,IB,DF)
C     --------------------------
C
      REAL D(*),DF(*)
      DATA SMALL /1.0E-20/
C
      IB=0
      IF (D(1).GE.0.0) IB=1
      IF (D(2).GE.0.0) IB=IB+2
      IF (D(3).GE.0.0) IB=IB+4
      IF (D(4).GE.0.0) IB=IB+8
      IF (IB.EQ.0 .OR. IB.EQ.15) RETURN
      DO 10 I=1,4
        J=1+MOD(I,4)
        IF (D(I)*D(J).LT.-SMALL) THEN
          DI=ABS(D(I))
          DF(I)=DI/(DI+ABS(D(J)))
        ENDIF
  10  CONTINUE
      END SUBROUTINE SBSRF1

      SUBROUTINE SBSRF2(XYZ,D1,D2,IB,FRC,VERT,EYE,LIGHT,LSHINE)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),D1(*),D2(*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      LOGICAL LSHINE
C
      IF (IB.EQ.15) THEN
        DO 10 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=XYZ(I)+D1(I)
          VERT(I,3)=VERT(I,2)+D2(I)
          VERT(I,4)=XYZ(I)+D2(I)
  10    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.1) THEN
        DO 20 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+FRC(1)*D1(I)
          VERT(I,3)=XYZ(I)+(1.0-FRC(4))*D2(I)
  20    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.2) THEN
        DO 30 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+FRC(2)*D2(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
  30    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.4) THEN
        DO 40 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=VERT(I,1)-FRC(3)*D1(I)
          VERT(I,3)=VERT(I,1)-(1.0-FRC(2))*D2(I)
  40    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.8) THEN
        DO 50 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=VERT(I,1)-FRC(4)*D2(I)
          VERT(I,3)=VERT(I,1)+(1.0-FRC(3))*D1(I)
  50    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.7) THEN
        DO 60 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+D1(I)
          VERT(I,3)=VERT(I,2)+D2(I)
          VERT(I,4)=VERT(I,3)-FRC(3)*D1(I)
          VERT(I,5)=XYZ(I)+(1.0-FRC(4))*D2(I)
  60    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.14) THEN
        DO 70 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+D2(I)
          VERT(I,3)=XYZ(I)+D2(I)
          VERT(I,4)=XYZ(I)+(1.0-FRC(4))*D2(I)
          VERT(I,5)=XYZ(I)+FRC(1)*D1(I)
  70    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.13) THEN
        DO 80 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=XYZ(I)+D2(I)
          VERT(I,3)=XYZ(I)
          VERT(I,4)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,5)=XYZ(I)+D1(I)+FRC(2)*D2(I)
  80    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.11) THEN
        DO 90 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=XYZ(I)
          VERT(I,3)=XYZ(I)+D1(I)
          VERT(I,4)=VERT(I,3)+FRC(2)*D2(I)
          VERT(I,5)=VERT(I,1)+(1.0-FRC(3))*D1(I)
  90    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.3) THEN
        DO 100 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=XYZ(I)+D1(I)
          VERT(I,3)=VERT(I,2)+FRC(2)*D2(I)
          VERT(I,4)=XYZ(I)+(1.0-FRC(4))*D2(I)
 100    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.6) THEN
        DO 110 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+D2(I)
          VERT(I,3)=VERT(I,2)-FRC(3)*D1(I)
          VERT(I,4)=XYZ(I)+FRC(1)*D1(I)
 110    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.12) THEN
        DO 120 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=XYZ(I)+D2(I)
          VERT(I,3)=VERT(I,2)-FRC(4)*D2(I)
          VERT(I,4)=XYZ(I)+D1(I)+FRC(2)*D2(I)
 120    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.9) THEN
        DO 130 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=XYZ(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,4)=VERT(I,1)+(1.0-FRC(3))*D1(I)
 130    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.5) THEN
        DO 140 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+FRC(1)*D1(I)
          VERT(I,3)=XYZ(I)+(1.0-FRC(4))*D2(I)
          VERT(I,4)=XYZ(I)+D1(I)+D2(I)
          VERT(I,5)=VERT(I,4)-FRC(3)*D1(I)
          VERT(I,6)=VERT(I,4)-(1.0-FRC(2))*D2(I)
 140    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
        CALL SBSRF6(EYE,3,VERT(1,4),LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.10) THEN
        DO 150 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+FRC(2)*D2(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,4)=XYZ(I)+D2(I)
          VERT(I,5)=VERT(I,4)-FRC(4)*D2(I)
          VERT(I,6)=VERT(I,4)+(1.0-FRC(3))*D1(I)
 150    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
        CALL SBSRF6(EYE,3,VERT(1,4),LSHINE,LIGHT,0) 
      ENDIF
      END SUBROUTINE SBSRF2

      SUBROUTINE SBSRF3(D,IVERT,DF,ISUMV,ISUMF)
C     -----------------------------------------
C
      REAL    D(*),DF(*)
      INTEGER IVERT(*),IC(8)
      DATA    SMALL /1.0E-20/
C
      ISUMV=0
      DO 10 I=1,8
        IF (D(I).LT.0.0) THEN
          IC(I)=0
        ELSE
          IC(I)=1
          ISUMV=ISUMV+1
        ENDIF
  10  CONTINUE
      IF (ISUMV.EQ.0 .OR. ISUMV.EQ.8) THEN
        ISUMV=0
        RETURN
      ENDIF
      IF (ISUMV.GT.4) THEN
        ISUMV=8-ISUMV
        DO 20 I=1,8
  20      IC(I)=MOD(IC(I)+1,2)
      ENDIF
      J=0
      DO 30 I=1,8
        IF (IC(I).EQ.1) THEN
          J=J+1
          IVERT(J)=I
        ENDIF
  30  CONTINUE
      ISUMF=0
      DO 40 I=1,4
        J=1+MOD(I,4)
        IF (D(I)*D(J).LT.-SMALL) THEN
          DI=ABS(D(I))
          DF(I)=DI/(DI+ABS(D(J)))
          ISUMF=ISUMF+1
        ENDIF
        K=I+4
        IF (D(I)*D(K).LT.-SMALL) THEN
          DI=ABS(D(I))
          DF(K)=DI/(DI+ABS(D(K)))
          ISUMF=ISUMF+1
        ENDIF
        L=J+4
        IF (D(K)*D(L).LT.-SMALL) THEN
          DK=ABS(D(K))
          DF(I+8)=DK/(DK+ABS(D(L)))
          ISUMF=ISUMF+1
        ENDIF
  40  CONTINUE
      END SUBROUTINE SBSRF3

      SUBROUTINE SBSRF4(DENS,N1,N2,N3,I0,J0,K0,GRDSCL,BAS,GRD)
C     --------------------------------------------------------
C
      REAL DENS(0:N1,0:N2,0:N3),GRDSCL(*),BAS(3,*),GRD(3,*)
      REAL G(3,8)
C
      IM=I0-1
      IF (IM.LT.0) IM=N1
      JM=J0-1
      IF (JM.LT.0) JM=N2
      KM=K0-1
      IF (KM.LT.0) KM=N3
      I1=I0+1
      J1=J0+1
      K1=K0+1
      IP=I1+1
      IF (IP.GT.N1) IP=0
      JP=J1+1
      IF (JP.GT.N2) JP=0
      KP=K1+1
      IF (KP.GT.N3) KP=0
      G(1,1)=GRDSCL(1)*(DENS(I1,J0,K0)-DENS(IM,J0,K0))
      G(2,1)=GRDSCL(2)*(DENS(I0,J1,K0)-DENS(I0,JM,K0))
      G(3,1)=GRDSCL(3)*(DENS(I0,J0,K1)-DENS(I0,J0,KM))
      G(1,2)=GRDSCL(1)*(DENS(IP,J0,K0)-DENS(I0,J0,K0))
      G(2,2)=GRDSCL(2)*(DENS(I1,J1,K0)-DENS(I1,JM,K0))
      G(3,2)=GRDSCL(3)*(DENS(I1,J0,K1)-DENS(I1,J0,KM))
      G(1,3)=GRDSCL(1)*(DENS(IP,J1,K0)-DENS(I0,J1,K0))
      G(2,3)=GRDSCL(2)*(DENS(I1,JP,K0)-DENS(I1,J0,K0))
      G(3,3)=GRDSCL(3)*(DENS(I1,J1,K1)-DENS(I1,J1,KM))
      G(1,4)=GRDSCL(1)*(DENS(I1,J1,K0)-DENS(IM,J1,K0))
      G(2,4)=GRDSCL(2)*(DENS(I0,JP,K0)-DENS(I0,J0,K0))
      G(3,4)=GRDSCL(3)*(DENS(I0,J1,K1)-DENS(I0,J1,KM))
      G(1,5)=GRDSCL(1)*(DENS(I1,J0,K1)-DENS(IM,J0,K1))
      G(2,5)=GRDSCL(2)*(DENS(I0,J1,K1)-DENS(I0,JM,K1))
      G(3,5)=GRDSCL(3)*(DENS(I0,J0,KP)-DENS(I0,J0,K0))
      G(1,6)=GRDSCL(1)*(DENS(IP,J0,K1)-DENS(I0,J0,K1))
      G(2,6)=GRDSCL(2)*(DENS(I1,J1,K1)-DENS(I1,JM,K1))
      G(3,6)=GRDSCL(3)*(DENS(I1,J0,KP)-DENS(I1,J0,K0))
      G(1,7)=GRDSCL(1)*(DENS(IP,J1,K1)-DENS(I0,J1,K1))
      G(2,7)=GRDSCL(2)*(DENS(I1,JP,K1)-DENS(I1,J0,K1))
      G(3,7)=GRDSCL(3)*(DENS(I1,J1,KP)-DENS(I1,J1,K0))
      G(1,8)=GRDSCL(1)*(DENS(I1,J1,K1)-DENS(IM,J1,K1))
      G(2,8)=GRDSCL(2)*(DENS(I0,JP,K1)-DENS(I0,J0,K1))
      G(3,8)=GRDSCL(3)*(DENS(I0,J1,KP)-DENS(I0,J1,K0))
      DO 20 J=1,8
        DO 10 I=1,3
  10      GRD(I,J)=G(1,J)*BAS(I,1)+G(2,J)*BAS(I,2)+G(3,J)*BAS(I,3)
  20  CONTINUE
      END SUBROUTINE SBSRF4

      SUBROUTINE SBSRF5(XYZ,DXYZ,ISV,ISF,IV,FRC,VERT,EYE,LIGHT,LSHINE)
C     ----------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV(*),IV4MAP(12)
      LOGICAL LSHINE
      DATA    IV4MAP /12,8,4,3,11,7,6,2,10,9,5,1/
C
      IF (ISV.EQ.1) THEN
        CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(1),EYE,LSHINE,LIGHT)
      ELSEIF (ISV.EQ.2) THEN
        IF (ISF.EQ.6) THEN
          CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(1),EYE,LSHINE,LIGHT)
          CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(2),EYE,LSHINE,LIGHT)
        ELSE
          IJDIF=IV(2)-IV(1)
          IF (IV(1).LE.4) THEN
            IF (IV(2).LE.4) THEN
              K2=IV(1)
              IF (IJDIF.EQ.3) K2=IV(2)
            ELSE
              K2=IV(2)
            ENDIF
          ELSE
            K2=IV(1)+4
            IF (IJDIF.EQ.3) K2=IV(2)+4
          ENDIF
          CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K2,EYE,LSHINE,LIGHT,1)
        ENDIF
      ELSEIF (ISV.EQ.3) THEN
        IF (ISF.EQ.9) THEN
          DO 10 I=1,3
  10        CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(I),EYE,LSHINE,LIGHT)
        ELSEIF (ISF.EQ.6) THEN
          DO 20 I1=1,3
            I2=1+MOD(I1,3)
            I=MIN(I1,I2)
            J=MAX(I1,I2)
            K2=0
            IJDIF=IV(J)-IV(I)
            IF (IV(I).LE.4) THEN
              IF (IV(J).LE.4) THEN
                IF (IJDIF.EQ.1) THEN
                  K2=IV(I)
                ELSEIF (IJDIF.EQ.3) THEN
                  K2=IV(J)
                ENDIF
              ELSE
                IF (IJDIF.EQ.4) K2=IV(J)
              ENDIF
            ELSE
              IF (IJDIF.EQ.1) THEN
                K2=IV(I)+4
              ELSEIF (IJDIF.EQ.3) THEN
                K2=IV(J)+4
              ENDIF
            ENDIF
            IF (K2.GT.0) GOTO 1
  20      CONTINUE
   1      CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K2,EYE,LSHINE,LIGHT,1)
        ELSE
          K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
          CALL SBSF5C(XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
        ENDIF
      ELSE
        IF (ISF.EQ.12) THEN
          DO 30 I=1,4
  30        CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(I),EYE,LSHINE,LIGHT)
        ELSEIF (ISF.EQ.4) THEN
          K4=(IV(1)+IV(2)+IV(3)+IV(4)-6)/4
          IF ((IV(2)-IV(1)).EQ.3) K4=6
          CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT,2)
        ELSEIF (ISF.EQ.6) THEN
          IF (IV(3).LE.4) THEN
            K3=IV(1)+IV(2)+IV(3)-6
            K4=MOD((IV(4)+K3),4)+3*K3
          ELSE
            IF (IV(2).GE.5) THEN
              K3=IV(2)+IV(3)+IV(4)-18
              K4=IV4MAP(MOD((IV(1)+K3),4)+3*K3)
            ELSE
              K4=12+IV(3)-IV(2)
              IF ((IV(1)+IV(2)+IV(3)+IV(4)).EQ.22) K4=29-K4
            ENDIF
          ENDIF
          CALL SBSF5D(XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT)
        ELSE
          K4=IV(1)+IV(2)+IV(3)+IV(4)
          IF (K4.EQ.16 .OR. K4.EQ.20) THEN          
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,IV(3),EYE,LSHINE,LIGHT,1)
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,IV(4),EYE,LSHINE,LIGHT,1)
          ELSEIF (K4.EQ.18) THEN
            K4A=IV(1)
            IF ((IV(2)-K4A).EQ.3) K4A=4
            K4B=9+MOD(K4A+1,4)
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K4A,EYE,LSHINE,LIGHT,1)
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K4B,EYE,LSHINE,LIGHT,1)
          ELSE
            IF (K4.EQ.14) THEN
              K4A=IV(4)
              K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
            ELSEIF (K4.EQ.22) THEN
              K4A=IV(1)
              K3=IV(2)+IV(3)+IV(4)-5+2*(IV(2)/5+2*(IV(3)/5)+IV(4)/5)
            ELSE
              IF (MOD((IV(1)+IV(2)),2).EQ.0) THEN
                IF (IV(4).EQ.6 .OR. (IV(3)-IV(2)).EQ.2) THEN
                  K4A=IV(2)
                  K3=IV(1)+IV(3)+IV(4)-5+2*(IV(1)/5+2*(IV(3)/5)+IV(4)/5)
                ELSE
                  K4A=IV(1)
                  K3=IV(2)+IV(3)+IV(4)-5+2*(IV(2)/5+2*(IV(3)/5)+IV(4)/5)
                ENDIF
              ELSE
                IF (IV(1).EQ.3 .OR. (IV(3)-IV(2)).EQ.2) THEN
                  K4A=IV(3)
                  K3=IV(1)+IV(2)+IV(4)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(4)/5)
                ELSE
                  K4A=IV(4)
                  K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
                ENDIF
              ENDIF
            ENDIF
            CALL SBSF5A(XYZ,DXYZ,FRC,VERT,K4A,EYE,LSHINE,LIGHT)
            CALL SBSF5C(XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
          ENDIF
        ENDIF
      ENDIF
      END SUBROUTINE SBSRF5

      SUBROUTINE SBSF5A(XYZ,DXYZ,FRC,VERT,IV,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      LOGICAL LSHINE
C
      IF (IV.LE.4) THEN
        J=IV
        K=1+MOD(IV+2,4)
        L=IV+4
      ELSE
        J=IV+4
        K=9+MOD(IV-2,4)
        L=IV
      ENDIF
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
  10  CONTINUE
      CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,1) 
      END SUBROUTINE SBSF5A

      SUBROUTINE SBSF5B(XYZ,DXYZ,FRC,VERT,KK,EYE,LSHINE,LIGHT,LL)
C     ------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IVL(4,12,2)
      LOGICAL LSHINE
      DATA    IVL /5,6,2,4,6,7,3,1,7,8,4,2,8,5,1,3,9,1,4,12,10,2,1,9,
     *     11,3,2,10,12,4,3,11,12,5,6,10,9,6,7,11,10,7,8,12,11,8,5,9,
     *     5,6,7,8,4,2,10,12,1,3,11,9,4,2,10,12,5,6,7,8,1,3,11,9,24*0/
C
      J=IVL(1,KK,LL)
      K=IVL(2,KK,LL)
      L=IVL(3,KK,LL)
      M=IVL(4,KK,LL)
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
        VERT(I,4)=XYZ(I)+DXYZ(I,M,1)+FRC(M)*DXYZ(I,M,2)
        VERT(I,5)=VERT(I,1)
        VERT(I,6)=0.25*(VERT(I,1)+VERT(I,2)+VERT(I,3)+VERT(I,4))
  10  CONTINUE
      DO 20 I=1,4
        CALL SBRCOP(VERT(1,I),VERT(1,7),6)
        CALL SBSRF6(EYE,3,VERT(1,6),LSHINE,LIGHT,1)
  20  CONTINUE
      END SUBROUTINE SBSF5B

      SUBROUTINE SBSF5C(XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV3(5,24)
      LOGICAL LSHINE
      DATA    IV3 /5,6,7,3,4,8,5,6,2,3,7,8,5,1,2,6,7,8,4,1,
     *             12,4,2,6,9,4,2,10,9,5,9,1,3,8,12,9,1,3,7,10,
     *             1,3,11,10,6,1,3,11,12,5,4,2,10,11,8,2,4,12,11,7,
     *             10,12,4,1,6,12,10,2,1,5,11,9,1,4,8,1,9,11,7,2,
     *             3,11,9,6,2,3,11,9,5,4,2,10,12,8,3,4,12,10,7,3,
     *             5,6,7,11,12,8,5,6,10,11,7,8,5,9,10,6,7,8,12,9/
C
      J=IV3(1,K3)
      K=IV3(2,K3)
      L=IV3(3,K3)
      M=IV3(4,K3)
      N=IV3(5,K3)
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
        VERT(I,4)=XYZ(I)+DXYZ(I,M,1)+FRC(M)*DXYZ(I,M,2)
        VERT(I,5)=XYZ(I)+DXYZ(I,N,1)+FRC(N)*DXYZ(I,N,2)
        VERT(I,6)=VERT(I,1)
        VERT(I,7)=0.2*(VERT(I,1)+VERT(I,2)+VERT(I,3)+VERT(I,4)+
     *                 VERT(I,5))
  10  CONTINUE
      DO 20 I=1,5
        CALL SBRCOP(VERT(1,I),VERT(1,8),6)
        CALL SBSRF6(EYE,3,VERT(1,7),LSHINE,LIGHT,1)
  20  CONTINUE
      END SUBROUTINE SBSF5C

      SUBROUTINE SBSF5D(XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV4(6,16)
      LOGICAL LSHINE
      DATA    IV4 /12,9,6,7,3,4,10,9,5,4,3,7,11,10,6,5,4,3,
     *             11,12,5,6,2,3,9,12,8,3,2,6,10,9,5,8,3,2,
     *             10,11,8,5,1,2,12,11,7,2,1,5,9,12,8,7,2,1,
     *             9,10,7,8,4,1,11,10,6,1,4,8,12,11,7,6,1,4,
     *             12,10,6,1,3,8,12,10,7,3,1,5,11,9,6,2,4,8,
     *             11,9,5,4,2,7/
      DATA    VNORM /0.1666666667/
C
      CALL SBRFIL(VERT(1,8),0.0,3)
      DO 20 J=1,6
        K=IV4(J,K4)
        DO 10 I=1,3
          VERT(I,J)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
          VERT(I,8)=VERT(I,8)+VERT(I,J)
  10    CONTINUE
  20  CONTINUE
      CALL SBRCOP(VERT,VERT(1,7),3)
      DO 30 I=1,3
  30    VERT(I,8)=VERT(I,8)*VNORM
      DO 40 I=1,6
        CALL SBRCOP(VERT(1,I),VERT(1,9),6)
        CALL SBSRF6(EYE,3,VERT(1,8),LSHINE,LIGHT,1)
  40  CONTINUE
      END SUBROUTINE SBSF5D

      SUBROUTINE SBSRF6(EYE,NV,VERT,LSHINE,LIGHT,INSIDE)
C     --------------------------------------------------
C
      REAL             EYE(*),VERT(3,*),LIGHT(*)
      LOGICAL          LSHINE
C
      REAL*8           XLNORM,ZZ,DZZ
      REAL             XW(20),YW(20),MTRX
      LOGICAL          LPS,LCOLOR
      COMMON  /SRFCOM/ GRDCUB(3,8),MTRX(3,3),ORIG(3),XL2,COL0,COLSCL
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C Carry out some initial checks and calculate the coordinates of the 
C projected triangle.
C
      IF (NV.LT.3 .OR. NV.GT.10) RETURN
      SMALL=1.0E-10
      XMIN=+1.0E20
      XMAX=-1.0E20
      YMIN=+1.0E20
      YMAX=-1.0E20
      DO 10 I=1,NV
        CALL SBLIN1(EYE,VERT(1,I),VERT(2,I),VERT(3,I),XW(I),YW(I))
        IF (XW(I).LT.XMIN) THEN
          XMIN=XW(I)
          ILEFT=I
        ENDIF
        IF (YW(I).LT.YMIN) THEN
          YMIN=YW(I)
          JBOTOM=I
        ENDIF
        XMAX=MAX(XW(I),XMAX)
        YMAX=MAX(YW(I),YMAX)
  10  CONTINUE
      IF (XMIN.GE.XTRC .OR. XMAX.LE.XBLC) RETURN
      IF (YMIN.GE.YTRC .OR. YMAX.LE.YBLC) RETURN
C
C Find the outward normal seen by the eye.
C
      AX=VERT(1,2)-VERT(1,1)
      AY=VERT(2,2)-VERT(2,1)
      AZ=VERT(3,2)-VERT(3,1)
      BX=VERT(1,1)-VERT(1,NV)
      BY=VERT(2,1)-VERT(2,NV)
      BZ=VERT(3,1)-VERT(3,NV)
      XN=BY*AZ-AY*BZ
      YN=BZ*AX-AZ*BX
      ZN=BX*AY-AX*BY
      TEN=XN*(EYE(1)-VERT(1,1))+YN*(EYE(2)-VERT(2,1))
     *   +ZN*(EYE(3)-VERT(3,1))
      IF (TEN.LT.0.0) THEN
        XN=-XN
        YN=-YN
        ZN=-ZN
        TEN=-TEN
      ENDIF
C
C Plot the projected triangle.
C
      XLNORM=DBLE(TEN)
      EYENRM=XN*EYE(1)+YN*EYE(2)+ZN*EYE(3)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      DYJ=1.0/DY
      DXI=1.0/DX
      SAFER=0.0001
      IF ((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
        JYMIN=INT((YMIN-YBLC)*DY)+2
        JYMAX=MIN(INT((YMAX-YBLC)*DY)+1,NYP)
        IF (JYMIN.GT.JYMAX) RETURN
        YJ=YBLC+(FLOAT(JYMIN-1)+SAFER)*DYJ
        NVL2=JBOTOM
        NVR2=JBOTOM
        J1=JYMIN
        DO 40 IVERT=1,NV
          IF (YJ.GT.YW(NVL2)) THEN
   1        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVL2)) GOTO 1
            YDIFL=YW(NVL2)-YW(NVL1)
            IF (ABS(YDIFL).LT.SMALL) YDIFL=SMALL
            GRADL=(XW(NVL2)-XW(NVL1))/YDIFL
          ENDIF
          IF (YJ.GT.YW(NVR2)) THEN
   2        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVR2)) GOTO 2
            YDIFR=YW(NVR2)-YW(NVR1)
            IF (ABS(YDIFR).LT.SMALL) YDIFR=SMALL
            GRADR=(XW(NVR2)-XW(NVR1))/YDIFR
          ENDIF
          IF (YW(NVL2).LT.YW(NVR2)) THEN
            J2=MIN(INT((YW(NVL2)-YBLC)*DY)+1,JYMAX)
          ELSE
            J2=MIN(INT((YW(NVR2)-YBLC)*DY)+1,JYMAX)
          ENDIF
          DO 30 J=J1,J2
            IF (J.GE.1) THEN
              XL=XW(NVL1)+GRADL*(YJ-YW(NVL1))
              XR=XW(NVR1)+GRADR*(YJ-YW(NVR1))
              ISTEP=1
              IX1=MAX(INT((XL-XBLC)*DX)+2,1)
              IX2=MIN(INT((XR-XBLC)*DX)+1,NXP)
              IF (IX1.GT.IX2) THEN
                ISTEP=-1
                IX1=MIN(IX1-1,NXP)
                IX2=MAX(IX2+1,1)
              ENDIF
              XI=XBLC+FLOAT(IX1-1)*DXI
              SDXI=FLOAT(ISTEP)*DXI
              DZZ=DBLE(SDXI*XN)
              ZZ=DBLE(EYENRM-XI*XN-YJ*YN)
              K=(J-1)*NXP+IX1
              DO 20 I=IX1,IX2,ISTEP
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))
                  IF (INSIDE.EQ.0) THEN
                    GX=XN
                    GY=YN
                    GZ=ZN
                  ELSE
                    CALL SBSF6A(X,Y,Z,ORIG,MTRX,GRDCUB,GX,GY,GZ)
                  ENDIF
                  CALL SBSF6B(EYE,X,Y,Z,GX,GY,GZ,LIGHT,XL2,LSHINE,CLR)
                  SBBUFF(KSTART+K)=COL0+COLSCL*CLR
                ENDIF
                XI=XI+SDXI
                ZZ=ZZ-DZZ
                K=K+ISTEP
  20          CONTINUE
            ENDIF
            YJ=YJ+DYJ
  30      CONTINUE
          J1=J2+1
          IF (J1.GT.JYMAX) RETURN
  40    CONTINUE
      ELSE
        IXMIN=INT((XMIN-XBLC)*DX)+2
        IXMAX=MIN(INT((XMAX-XBLC)*DX)+1,NXP)
        IF (IXMIN.GT.IXMAX) RETURN
        XI=XBLC+(FLOAT(IXMIN-1)+SAFER)*DXI
        NVL2=ILEFT
        NVR2=ILEFT
        I1=IXMIN
        DO 70 IVERT=1,NV
          IF (XI.GT.XW(NVL2)) THEN
   3        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVL2)) GOTO 3
            XDIFL=XW(NVL2)-XW(NVL1)
            IF (ABS(XDIFL).LT.SMALL) XDIFL=SMALL
            GRADL=(YW(NVL2)-YW(NVL1))/XDIFL
          ENDIF
          IF (XI.GT.XW(NVR2)) THEN
   4        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVR2)) GOTO 4
            XDIFR=XW(NVR2)-XW(NVR1)
            IF (ABS(XDIFR).LT.SMALL) XDIFR=SMALL
            GRADR=(YW(NVR2)-YW(NVR1))/XDIFR
          ENDIF
          IF (XW(NVL2).LT.XW(NVR2)) THEN
            I2=MIN(INT((XW(NVL2)-XBLC)*DX)+1,IXMAX)
          ELSE
            I2=MIN(INT((XW(NVR2)-XBLC)*DX)+1,IXMAX)
          ENDIF
          DO 60 I=I1,I2
            IF (I.GE.1) THEN
              YL=YW(NVL1)+GRADL*(XI-XW(NVL1))
              YR=YW(NVR1)+GRADR*(XI-XW(NVR1))
              ISTEP=1
              JY1=MAX(INT((YL-YBLC)*DY)+2,1)
              JY2=MIN(INT((YR-YBLC)*DY)+1,NYP)
              IF (JY1.GT.JY2) THEN
                ISTEP=-1
                JY1=MIN(JY1-1,NYP)
                JY2=MAX(JY2+1,1)
              ENDIF
              YJ=YBLC+FLOAT(JY1-1)*DYJ
              SDYJ=FLOAT(ISTEP)*DYJ
              DZZ=DBLE(SDYJ*YN)
              ZZ=DBLE(EYENRM-YJ*YN-XI*XN)
              K=(JY1-1)*NXP+I
              KSTEP=ISTEP*NXP
              DO 50 J=JY1,JY2,ISTEP
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))
                  IF (INSIDE.EQ.0) THEN
                    GX=XN
                    GY=YN
                    GZ=ZN
                  ELSE
                    CALL SBSF6A(X,Y,Z,ORIG,MTRX,GRDCUB,GX,GY,GZ)
                  ENDIF
                  CALL SBSF6B(EYE,X,Y,Z,GX,GY,GZ,LIGHT,XL2,LSHINE,CLR)
                  SBBUFF(KSTART+K)=COL0+COLSCL*CLR
                ENDIF
                YJ=YJ+SDYJ
                ZZ=ZZ-DZZ
                K=K+KSTEP
  50          CONTINUE
            ENDIF
            XI=XI+DXI
  60      CONTINUE
          I1=I2+1
          IF (I1.GT.IXMAX) RETURN
  70    CONTINUE
      ENDIF
      END SUBROUTINE SBSRF6

      SUBROUTINE SBSF6A(X,Y,Z,ORIG,MTRX,GRD,XN,YN,ZN)
C     -----------------------------------------------
C
      REAL ORIG(*),MTRX(3,*),GRD(3,*)
      DATA ZERO,ONE /0.00001,0.99999/
C
      X0=X-ORIG(1)
      Y0=Y-ORIG(2)
      Z0=Z-ORIG(3)
      XI=X0*MTRX(1,1)+Y0*MTRX(2,1)+Z0*MTRX(3,1)
      YJ=X0*MTRX(1,2)+Y0*MTRX(2,2)+Z0*MTRX(3,2)
      ZK=X0*MTRX(1,3)+Y0*MTRX(2,3)+Z0*MTRX(3,3)
      DX=MIN(ONE,MAX(ZERO,XI-FLOAT(INT(XI))))
      DY=MIN(ONE,MAX(ZERO,YJ-FLOAT(INT(YJ))))
      DZ=MIN(ONE,MAX(ZERO,ZK-FLOAT(INT(ZK))))
      XN1=(1.0-DX)*(GRD(1,1)+DY*(GRD(1,4)-GRD(1,1)))
     *        +DX*(GRD(1,2)+DY*(GRD(1,3)-GRD(1,2)))
      XN2=(1.0-DX)*(GRD(1,5)+DY*(GRD(1,8)-GRD(1,5)))
     *        +DX*(GRD(1,6)+DY*(GRD(1,7)-GRD(1,6)))
      XN=XN1+DZ*(XN2-XN1)
      YN1=(1.0-DX)*(GRD(2,1)+DY*(GRD(2,4)-GRD(2,1)))
     *        +DX*(GRD(2,2)+DY*(GRD(2,3)-GRD(2,2)))
      YN2=(1.0-DX)*(GRD(2,5)+DY*(GRD(2,8)-GRD(2,5)))
     *        +DX*(GRD(2,6)+DY*(GRD(2,7)-GRD(2,6)))
      YN=YN1+DZ*(YN2-YN1)
      ZN1=(1.0-DX)*(GRD(3,1)+DY*(GRD(3,4)-GRD(3,1)))
     *        +DX*(GRD(3,2)+DY*(GRD(3,3)-GRD(3,2)))
      ZN2=(1.0-DX)*(GRD(3,5)+DY*(GRD(3,8)-GRD(3,5)))
     *        +DX*(GRD(3,6)+DY*(GRD(3,7)-GRD(3,6)))
      ZN=ZN1+DZ*(ZN2-ZN1)
      END SUBROUTINE SBSF6A

      SUBROUTINE SBSF6B(EYE,X,Y,Z,XN,YN,ZN,LIGHT,XL2,LSHINE,COLOUR)
C     -------------------------------------------------------------
C
      REAL    EYE(*),LIGHT(*)
      LOGICAL LSHINE
      DATA    SMALL2 /1.0E-20/
C
      COLOUR=0.0
      XNL=XN*LIGHT(1)+YN*LIGHT(2)+ZN*LIGHT(3)
      IF (XNL.GE.0.0) RETURN
      XN2=XN**2+YN**2+ZN**2+SMALL2
      IF (LSHINE) THEN
        RFNORM=2.0*XNL/XN2
        RX=LIGHT(1)-XN*RFNORM
        RY=LIGHT(2)-YN*RFNORM
        RZ=LIGHT(3)-ZN*RFNORM
        VX=EYE(1)-X
        VY=EYE(2)-Y
        VZ=EYE(3)-Z
        XRV=RX*VX+RY*VY+RZ*VZ
        IF (XRV.LT.0.0) RETURN
        V2=VX**2+VY**2+VZ**2
        COLOUR=MIN(XRV**2/(ABS(XL2*V2)+SMALL2),1.0)
      ELSE
        COLOUR=MIN(-XNL/SQRT(ABS(XL2*XN2)+SMALL2),1.0)
      ENDIF
      END SUBROUTINE SBSF6B

      SUBROUTINE SB2SRF(EYE,LATICE,DENS,N1,N2,DLOW,DHIGH,DVERT,IC1,IC2,
     *                  NCBAND,LIGHT,LSHINE)
C     -----------------------------------------------------------------
C
      REAL             EYE(*),LATICE(3,*),DENS(0:N1,0:N2),LIGHT(*)
      LOGICAL          LSHINE
C
      REAL*8           XLNORM,ZZ,DZZ,DSMAL2
      REAL             LATIC2(3,4),BAS(3,3),MTRX(3,3)
      REAL             VERT(3,8),XW(20),YW(20),GRDSCL(3)
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots a 3-d surface given a 2-d unit-cell
C    of density. All (x,y,z) values are taken to be given in world 
C    coordinates. The z-component of the eye-poisition should be 
C    positive and that of all the lattice-vertices should be negative; 
C    the viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    LATICE   R*4    I     3 x 3    (x,y,z) coordinates of the origin
C                                   and the a and b lattice-vertices.
C    DENS     R*4    I     (N1+1)   The density at regular points within
C                        x (N2+1)   the unit cell, wrapped around so
C                                   that DENS(0,J)=DENS(N1,J) etc..
C    N1,N2    I*4    I       -      The dimensions of the unit-cell grid.
C    DLOW     R*4    I       -      Lowest density to be plotted.
C    DHIGH    R*4    I       -      Highest density to be plotted.
C    DVERT    R*4    I       -      "Vertical" world-coordinate length
C                                   corresponding to density-range.
C    IC1,IC2  I*4    I       -      Lowest and highest colour-index to
C                                   be used for the rendering.
C    NCBAND   I*4    I       -      Number of colour-bands for the
C                                   height, so that the number of shades
C                                   per band = (IC2-IC1+1)/NCBAND.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    LSHINE   L*1    I       -      Shiny surface if TRUE, else diffuse.
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBLIN1     Calculates the projection of (x,y,z) on viewing screen.
C     SB2SR1     Calculates vertecies of triangular breakdown of grid.
C     SB2SR3     Calculates the normal of the surface.
C     SB2SR4     Calculates the appropriate colour for a given pixel.
C
C History
C   D. S. Sivia        1 Jun 1995  Initial release.
C   D. S. Sivia        7 Jul 1995  Fixed bug in determinant calculation.
C   D. S. Sivia       20 Oct 1995  Speeded up computations slightly.
C   D. S. Sivia       26 Oct 1995  Completely new algorithm!
C   D. S. Sivia       24 Oct 1997 "Safe-guarded" some rounding errors.
C-----------------------------------------------------------------------
C
C Carry out some initial checks.
C
      BIG=1.0E20
      SMALL=1.0E-10
      SMALL2=SMALL**2
      DSMAL2=DBLE(SMALL2)
      IF (EYE(3).LE.SMALL) RETURN
      IF (N1.LT.1 .OR. N2.LT.1) RETURN
      DRANGE=DHIGH-DLOW
      IF (DRANGE.LE.SMALL) RETURN
      XL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
      IF (XL2.LT.SMALL) RETURN
      DO 10 I=1,3
        BAS(I,1)=LATICE(I,2)-LATICE(I,1)
        BAS(I,2)=LATICE(I,3)-LATICE(I,1)
  10  CONTINUE
      CX=BAS(2,1)*BAS(3,2)-BAS(2,2)*BAS(3,1)
      CY=BAS(3,1)*BAS(1,2)-BAS(3,2)*BAS(1,1)
      CZ=BAS(1,1)*BAS(2,2)-BAS(1,2)*BAS(2,1)
      CSCL=DVERT/SQRT(CX**2+CY**2+CZ**2+SMALL**2)
      BAS(1,3)=CX*CSCL
      BAS(2,3)=CY*CSCL
      BAS(3,3)=CZ*CSCL
      DMIN=MIN(DLOW,0.0)
      DMAX=MAX(DHIGH,0.0)
      ORGSCL=DMIN/(DMAX-DMIN+SMALL)
      DO 11 I=1,3
        LATIC2(I,1)=LATICE(I,1)+ORGSCL*BAS(I,3)
        LATIC2(I,2)=LATIC2(I,1)+BAS(I,1)
        LATIC2(I,3)=LATIC2(I,1)+BAS(I,2)
        LATIC2(I,4)=LATIC2(I,1)+BAS(I,3)
  11  CONTINUE
      IF (LATIC2(3,1).GE.0.0) RETURN
      ZFAR=LATIC2(3,2)+LATIC2(3,3)+LATIC2(3,4)-2.0*LATIC2(3,1)
      IF (ZFAR.GE.0.0) RETURN
      DO 12 J=1,3
        IF (LATIC2(3,J+1).GE.0.0) RETURN
        IF ((ZFAR-BAS(3,J)).GE.0.0) RETURN
        BAS2J=BAS(1,J)**2+BAS(2,J)**2+BAS(3,J)**2
        IF (BAS2J.LT.SMALL2) RETURN
  12  CONTINUE
      EYLAT1=EYE(1)-LATIC2(1,1)
      EYLAT2=EYE(2)-LATIC2(2,1)
      NTOT=(N1+1)*(N2+1)
C
C Set up matrix for real-space to lattice-index transformation.
C
      SAFE=0.0001
      XN1=0.99999*FLOAT(N1)
      XN2=0.99999*FLOAT(N2)
      XN3=0.99999
      DET=BAS(1,1)*BAS(2,2)*BAS(3,3)+BAS(1,2)*BAS(2,3)*BAS(3,1)
     *   +BAS(1,3)*BAS(2,1)*BAS(3,2)-BAS(3,1)*BAS(2,2)*BAS(1,3)
     *   -BAS(3,2)*BAS(2,3)*BAS(1,1)-BAS(3,3)*BAS(2,1)*BAS(1,2)
      IF (ABS(DET).LT.SMALL2) RETURN
      DETNRM=1.0/DET
      MTRX(1,1)=XN1*DETNRM*(BAS(2,2)*BAS(3,3)-BAS(2,3)*BAS(3,2))
      MTRX(1,2)=XN2*DETNRM*(BAS(2,3)*BAS(3,1)-BAS(2,1)*BAS(3,3))
      MTRX(1,3)=XN3*DETNRM*(BAS(2,1)*BAS(3,2)-BAS(2,2)*BAS(3,1))
      MTRX(2,1)=XN1*DETNRM*(BAS(3,2)*BAS(1,3)-BAS(3,3)*BAS(1,2))
      MTRX(2,2)=XN2*DETNRM*(BAS(3,3)*BAS(1,1)-BAS(3,1)*BAS(1,3))
      MTRX(2,3)=XN3*DETNRM*(BAS(3,1)*BAS(1,2)-BAS(3,2)*BAS(1,1))
      MTRX(3,1)=XN1*DETNRM*(BAS(1,2)*BAS(2,3)-BAS(1,3)*BAS(2,2))
      MTRX(3,2)=XN2*DETNRM*(BAS(1,3)*BAS(2,1)-BAS(1,1)*BAS(2,3))
      MTRX(3,3)=XN3*DETNRM*(BAS(1,1)*BAS(2,2)-BAS(1,2)*BAS(2,1))
C
C Some general initialisations.
C
      GRDSCL(1)=-FLOAT(N1)
      GRDSCL(2)=-FLOAT(N2)
      GRDSCL(3)=1.0/DRANGE
      IF (IC2.LT.IC1) RETURN
      NSHADS=MAX((IC2-IC1+1)/MAX(NCBAND,1),1)
      COLSCL=FLOAT(NSHADS-1)
      COL0=FLOAT(IC1)
      DKSCL=0.9999*FLOAT((IC2-IC1+1)/NSHADS)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      DYJ=1.0/DY
      DXI=1.0/DX
C
C New algorithm here: divide each box of 4 grid-points into four 
C triangles and paint them.
C
      DO 70 JLAT=1,N2
        DO 60 ILAT=1,N1
          CALL SB2SR1(ILAT,JLAT,DENS,N1,N2,DLOW,DRANGE,BAS,LATIC2,VERT)
          DO 50 IV=1,4
            CALL SBRCOP(VERT(1,IV),VERT(1,6),6)
            XMIN=+BIG
            XMAX=-BIG
            YMIN=+BIG
            YMAX=-BIG
            DO 20 I=1,3
              II=5+I
              CALL SBLIN1(EYE,VERT(1,II),VERT(2,II),VERT(3,II),
     *                    XW(I),YW(I))
              IF (XW(I).LT.XMIN) THEN
                XMIN=XW(I)
                ILEFT=I
              ENDIF
              IF (YW(I).LT.YMIN) THEN
                YMIN=YW(I)
                JBOTOM=I
              ENDIF
              XMAX=MAX(XW(I),XMAX)
              YMAX=MAX(YW(I),YMAX)
  20        CONTINUE
            IF (XMIN.GE.XTRC .OR. XMAX.LE.XBLC) GOTO 50
            IF (YMIN.GE.YTRC .OR. YMAX.LE.YBLC) GOTO 50
            AX=VERT(1,7)-VERT(1,6)
            AY=VERT(2,7)-VERT(2,6)
            AZ=VERT(3,7)-VERT(3,6)
            BX=VERT(1,6)-VERT(1,8)
            BY=VERT(2,6)-VERT(2,8)
            BZ=VERT(3,6)-VERT(3,8)
            XN=BY*AZ-AY*BZ
            YN=BZ*AX-AZ*BX
            ZN=BX*AY-AX*BY
            TEN=XN*(EYE(1)-VERT(1,6))+YN*(EYE(2)-VERT(2,6))
     *         +ZN*(EYE(3)-VERT(3,6))
            XLNORM=DBLE(TEN)
            EYENRM=XN*EYE(1)+YN*EYE(2)+ZN*EYE(3)
            SAFER=0.0001
            IF ((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
              JYMIN=INT((YMIN-YBLC)*DY)+2
              JYMAX=MIN(INT((YMAX-YBLC)*DY)+1,NYP)
              IF (JYMIN.GT.JYMAX) GOTO 50
              YJ=YBLC+(FLOAT(JYMIN-1)+SAFER)*DYJ
              NVL2=JBOTOM
              NVR2=JBOTOM
              J1=JYMIN
              DO 30 IVERT=1,3
                IF (YJ.GT.YW(NVL2)) THEN
   1              NVL1=NVL2
                  NVL2=NVL1-1
                  IF (NVL2.LT.1) NVL2=3
                  IF (NVL2.EQ.JBOTOM) GOTO 50
                  IF (YJ.GT.YW(NVL2)) GOTO 1
                  YDIFL=YW(NVL2)-YW(NVL1)
                  IF (ABS(YDIFL).LT.SMALL) YDIFL=SMALL
                  GRADL=(XW(NVL2)-XW(NVL1))/YDIFL
                ENDIF
                IF (YJ.GT.YW(NVR2)) THEN
   2              NVR1=NVR2
                  NVR2=NVR1+1
                  IF (NVR2.GT.3) NVR2=1
                  IF (NVR2.EQ.JBOTOM) GOTO 50
                  IF (YJ.GT.YW(NVR2)) GOTO 2
                  YDIFR=YW(NVR2)-YW(NVR1)
                  IF (ABS(YDIFR).LT.SMALL) YDIFR=SMALL
                  GRADR=(XW(NVR2)-XW(NVR1))/YDIFR
                ENDIF
                IF (YW(NVL2).LT.YW(NVR2)) THEN
                  J2=MIN(INT((YW(NVL2)-YBLC)*DY)+1,JYMAX)
                ELSE
                  J2=MIN(INT((YW(NVR2)-YBLC)*DY)+1,JYMAX)
                ENDIF
                DO 29 J=J1,J2
                  IF (J.GE.1) THEN
                    XL=XW(NVL1)+GRADL*(YJ-YW(NVL1))
                    XR=XW(NVR1)+GRADR*(YJ-YW(NVR1))
                    ISTEP=1
                    IX1=MAX(INT((XL-XBLC)*DX)+2,1)
                    IX2=MIN(INT((XR-XBLC)*DX)+1,NXP)
                    IF (IX1.GT.IX2) THEN
                      ISTEP=-1
                      IX1=MIN(IX1-1,NXP)
                      IX2=MAX(IX2+1,1)
                    ENDIF
                    XI=XBLC+FLOAT(IX1-1)*DXI
                    DXISTP=FLOAT(ISTEP)*DXI
                    DZZ=DBLE(DXISTP*XN)
                    ZZ=DBLE(EYENRM-XI*XN-YJ*YN)
                    K=(J-1)*NXP+IX1
                    DO 28 I=IX1,IX2,ISTEP
                      XLAMDA=SNGL(XLNORM/ZZ)
                      Z=EYE(3)*(1.0-XLAMDA)
                      IF (Z.GT.SBBUFF(K)) THEN
                        SBBUFF(K)=Z
                        IF (IC1.EQ.IC2) THEN
                          SBBUFF(K+KSTART)=COL0
                        ELSE
                          XDI=EYLAT1+XLAMDA*(XI-EYE(1))
                          YDJ=EYLAT2+XLAMDA*(YJ-EYE(2))
                          ZDK=Z-LATIC2(3,1)
                          IF (NCBAND.GT.1) THEN
                            DK=MAX(MIN((XDI*MTRX(1,3)+YDJ*MTRX(2,3)
     *                        +ZDK*MTRX(3,3)),XN3),SAFE)
                            COL0=FLOAT(IC1+NSHADS*INT(DK*DKSCL))
                          ENDIF
                          DI=MAX(MIN((XDI*MTRX(1,1)+YDJ*MTRX(2,1)
     *                      +ZDK*MTRX(3,1)),XN1),SAFE)
                          DJ=MAX(MIN((XDI*MTRX(1,2)+YDJ*MTRX(2,2)
     *                      +ZDK*MTRX(3,2)),XN2),SAFE)
                          CALL SB2SR3(DENS,N1,N2,DI,DJ,BAS,GX,GY,GZ,
     *                                GRDSCL)
                          CALL SB2SR4(EYE,XI,YJ,0.0,GX,GY,GZ,LIGHT,
     *                                XL2,SMALL2,LSHINE,COLOUR)
                          SBBUFF(K+KSTART)=COL0+COLOUR*COLSCL
                        ENDIF
                      ENDIF
                      XI=XI+DXISTP
                      ZZ=ZZ-DZZ
                      K=K+ISTEP
  28                CONTINUE
                  ENDIF
                  YJ=YJ+DYJ
  29            CONTINUE
                J1=J2+1
                IF (J1.GT.JYMAX) GOTO 50
  30          CONTINUE
            ELSE
              IXMIN=INT((XMIN-XBLC)*DX)+2
              IXMAX=MIN(INT((XMAX-XBLC)*DX)+1,NXP)
              IF (IXMIN.GT.IXMAX) GOTO 50
              XI=XBLC+(FLOAT(IXMIN-1)+SAFER)*DXI
              NVL2=ILEFT
              NVR2=ILEFT
              I1=IXMIN
              DO 40 IVERT=1,3
                IF (XI.GT.XW(NVL2)) THEN
   3              NVL1=NVL2
                  NVL2=NVL1-1
                  IF (NVL2.LT.1) NVL2=3
                  IF (NVL2.EQ.ILEFT) GOTO 50
                  IF (XI.GT.XW(NVL2)) GOTO 3
                  XDIFL=XW(NVL2)-XW(NVL1)
                  IF (ABS(XDIFL).LT.SMALL) XDIFL=SMALL
                  GRADL=(YW(NVL2)-YW(NVL1))/XDIFL
                ENDIF
                IF (XI.GT.XW(NVR2)) THEN
   4              NVR1=NVR2
                  NVR2=NVR1+1
                  IF (NVR2.GT.3) NVR2=1
                  IF (NVR2.EQ.ILEFT) GOTO 50
                  IF (XI.GT.XW(NVR2)) GOTO 4
                  XDIFR=XW(NVR2)-XW(NVR1)
                  IF (ABS(XDIFR).LT.SMALL) XDIFR=SMALL
                  GRADR=(YW(NVR2)-YW(NVR1))/XDIFR
               ENDIF
                IF (XW(NVL2).LT.XW(NVR2)) THEN
                  I2=MIN(INT((XW(NVL2)-XBLC)*DX)+1,IXMAX)
                ELSE
                  I2=MIN(INT((XW(NVR2)-XBLC)*DX)+1,IXMAX)
                ENDIF
                DO 39 I=I1,I2
                  IF (I.GE.1) THEN
                    YL=YW(NVL1)+GRADL*(XI-XW(NVL1))
                    YR=YW(NVR1)+GRADR*(XI-XW(NVR1))
                    JSTEP=1
                    JY1=MAX(INT((YL-YBLC)*DY)+2,1)
                    JY2=MIN(INT((YR-YBLC)*DY)+1,NYP)
                    IF (JY1.GT.JY2) THEN
                      JSTEP=-1
                      JY1=MIN(JY1-1,NYP)
                      JY2=MAX(JY2+1,1)
                    ENDIF
                    YJ=YBLC+FLOAT(JY1-1)*DYJ
                    DYJSTP=FLOAT(JSTEP)*DYJ
                    DZZ=DBLE(DYJSTP*YN)
                    ZZ=DBLE(EYENRM-YJ*YN-XI*XN)
                    K=(JY1-1)*NXP+I
                    KSTEP=JSTEP*NXP
                    DO 38 J=JY1,JY2,JSTEP
                      XLAMDA=SNGL(XLNORM/ZZ)
                      Z=EYE(3)*(1.0-XLAMDA)
                      IF (Z.GT.SBBUFF(K)) THEN
                        SBBUFF(K)=Z
                        IF (IC1.EQ.IC2) THEN
                          SBBUFF(K+KSTART)=COL0
                        ELSE
                          XDI=EYLAT1+XLAMDA*(XI-EYE(1))
                          YDJ=EYLAT2+XLAMDA*(YJ-EYE(2))
                          ZDK=Z-LATIC2(3,1)
                          IF (NCBAND.GT.1) THEN
                            DK=MAX(MIN((XDI*MTRX(1,3)+YDJ*MTRX(2,3)
     *                        +ZDK*MTRX(3,3)),XN3),SAFE)
                            COL0=FLOAT(IC1+NSHADS*INT(DK*DKSCL))
                          ENDIF
                          DI=MAX(MIN((XDI*MTRX(1,1)+YDJ*MTRX(2,1)
     *                      +ZDK*MTRX(3,1)),XN1),SAFE)
                          DJ=MAX(MIN((XDI*MTRX(1,2)+YDJ*MTRX(2,2)
     *                      +ZDK*MTRX(3,2)),XN2),SAFE)
                          CALL SB2SR3(DENS,N1,N2,DI,DJ,BAS,GX,GY,GZ,
     *                                GRDSCL)
                          CALL SB2SR4(EYE,XI,YJ,0.0,GX,GY,GZ,LIGHT,
     *                                XL2,SMALL2,LSHINE,COLOUR)
                          SBBUFF(K+KSTART)=COL0+COLOUR*COLSCL
                        ENDIF
                      ENDIF
                      YJ=YJ+DYJSTP
                      ZZ=ZZ-DZZ
                      K=K+KSTEP
  38                CONTINUE
                  ENDIF
                  XI=XI+DXI
  39            CONTINUE
                I1=I2+1
                IF (I1.GT.IXMAX) GOTO 50
  40          CONTINUE
            ENDIF
  50      CONTINUE
  60    CONTINUE
  70  CONTINUE
      END       SUBROUTINE SB2SRF

      SUBROUTINE SB2SR1(I1,J1,DENS,N1,N2,DLOW,DRANGE,BAS,LATICE,VERT)
C     ---------------------------------------------------------------
C
      REAL DENS(0:N1,0:N2),BAS(3,*),LATICE(*),VERT(3,*)
C
      I0=I1-1
      J0=J1-1
      XNORM=1.0/FLOAT(N1)
      YNORM=1.0/FLOAT(N2)
      ZNORM=1.0/DRANGE
      D00=MAX(MIN((DENS(I0,J0)-DLOW)*ZNORM,1.0),0.0)
      D10=MAX(MIN((DENS(I1,J0)-DLOW)*ZNORM,1.0),0.0)
      D11=MAX(MIN((DENS(I1,J1)-DLOW)*ZNORM,1.0),0.0)
      D01=MAX(MIN((DENS(I0,J1)-DLOW)*ZNORM,1.0),0.0)
      X0=FLOAT(I0)*XNORM
      X1=FLOAT(I1)*XNORM
      Y0=FLOAT(J0)*YNORM
      Y1=FLOAT(J1)*YNORM
      VERT(1,1)=X0*BAS(1,1)+Y0*BAS(1,2)+D00*BAS(1,3)+LATICE(1)
      VERT(2,1)=X0*BAS(2,1)+Y0*BAS(2,2)+D00*BAS(2,3)+LATICE(2)
      VERT(3,1)=X0*BAS(3,1)+Y0*BAS(3,2)+D00*BAS(3,3)+LATICE(3)
      VERT(1,2)=X1*BAS(1,1)+Y0*BAS(1,2)+D10*BAS(1,3)+LATICE(1)
      VERT(2,2)=X1*BAS(2,1)+Y0*BAS(2,2)+D10*BAS(2,3)+LATICE(2)
      VERT(3,2)=X1*BAS(3,1)+Y0*BAS(3,2)+D10*BAS(3,3)+LATICE(3)
      VERT(1,3)=X1*BAS(1,1)+Y1*BAS(1,2)+D11*BAS(1,3)+LATICE(1)
      VERT(2,3)=X1*BAS(2,1)+Y1*BAS(2,2)+D11*BAS(2,3)+LATICE(2)
      VERT(3,3)=X1*BAS(3,1)+Y1*BAS(3,2)+D11*BAS(3,3)+LATICE(3)
      VERT(1,4)=X0*BAS(1,1)+Y1*BAS(1,2)+D01*BAS(1,3)+LATICE(1)
      VERT(2,4)=X0*BAS(2,1)+Y1*BAS(2,2)+D01*BAS(2,3)+LATICE(2)
      VERT(3,4)=X0*BAS(3,1)+Y1*BAS(3,2)+D01*BAS(3,3)+LATICE(3)
      VERT(1,5)=VERT(1,1)
      VERT(2,5)=VERT(2,1)
      VERT(3,5)=VERT(3,1)
      XM=0.5*FLOAT(I0+I1)*XNORM
      YM=0.5*FLOAT(J0+J1)*YNORM
      DM=0.25*(DENS(I0,J0)+DENS(I1,J0)+DENS(I1,J1)+DENS(I0,J1))
      DM=MAX(MIN((DM-DLOW)*ZNORM,1.0),0.0)
      VERT(1,8)=XM*BAS(1,1)+YM*BAS(1,2)+DM*BAS(1,3)+LATICE(1)
      VERT(2,8)=XM*BAS(2,1)+YM*BAS(2,2)+DM*BAS(2,3)+LATICE(2)
      VERT(3,8)=XM*BAS(3,1)+YM*BAS(3,2)+DM*BAS(3,3)+LATICE(3)
      END SUBROUTINE SB2SR1

      SUBROUTINE SB2SR3(DENS,N1,N2,XI,YJ,BAS,XN,YN,ZN,GRDSCL)
C     -------------------------------------------------------
C
      REAL DENS(0:N1,0:N2),BAS(3,*),GRDSCL(*)
C
      I=INT(XI)
      J=INT(YJ)
      II=I+1
      JJ=J+1
      DX=XI-FLOAT(I)
      DY=YJ-FLOAT(J)
      CALL S2SR3A(I,DX,N1,IM,I0,IP,DDX)
      XM=DENS(IM,J)+DY*(DENS(IM,JJ)-DENS(IM,J))
      X0=DENS(I0,J)+DY*(DENS(I0,JJ)-DENS(I0,J))
      XP=DENS(IP,J)+DY*(DENS(IP,JJ)-DENS(IP,J))
      GX=GRDSCL(1)*GRDSCL(3)*((1.0-DDX)*(X0-XM)+DDX*(XP-X0))
      CALL S2SR3A(J,DY,N2,JM,J0,JP,DDY)
      YM=DENS(I,JM)+DX*(DENS(II,JM)-DENS(I,JM))
      Y0=DENS(I,J0)+DX*(DENS(II,J0)-DENS(I,J0))
      YP=DENS(I,JP)+DX*(DENS(II,JP)-DENS(I,JP))
      GY=GRDSCL(2)*GRDSCL(3)*((1.0-DDY)*(Y0-YM)+DDY*(YP-Y0))
      GZ=1.0
      XN=GX*BAS(1,1)+GY*BAS(1,2)+GZ*BAS(1,3)
      YN=GX*BAS(2,1)+GY*BAS(2,2)+GZ*BAS(2,3)
      ZN=GX*BAS(3,1)+GY*BAS(3,2)+GZ*BAS(3,3)
      END SUBROUTINE SB2SR3

      SUBROUTINE S2SR3A(I,DX,NX,IM,I0,IP,DDX)
C     ---------------------------------------
C
      IF (DX.LT.0.5) THEN
        DDX=DX+0.5
        I0=I
        IM=I0-1
        IP=I0+1
        IF (IM.LT.0) IM=NX
      ELSE
        DDX=DX-0.5
        I0=I+1
        IM=I0-1
        IP=I0+1
        IF (IP.GT.NX) IP=0
      ENDIF
      END SUBROUTINE S2SR3A

      SUBROUTINE SB2SR4(EYE,X,Y,Z,XN,YN,ZN,LIGHT,XL2,SMALL2,LSHINE,
     *                  COLOUR)
C     -------------------------------------------------------------
C
      REAL    EYE(*),LIGHT(*)
      LOGICAL LSHINE
C
      COLOUR=0.0
      XNL=XN*LIGHT(1)+YN*LIGHT(2)+ZN*LIGHT(3)
      IF (XNL.GE.0.0) RETURN
      XN2=XN**2+YN**2+ZN**2+SMALL2
      IF (LSHINE) THEN
        RFNORM=2.0*XNL/XN2
        RX=LIGHT(1)-XN*RFNORM
        RY=LIGHT(2)-YN*RFNORM
        RZ=LIGHT(3)-ZN*RFNORM
        VX=EYE(1)-X
        VY=EYE(2)-Y
        VZ=EYE(3)-Z
        XRV=RX*VX+RY*VY+RZ*VZ
        IF (XRV.LT.0.0) RETURN
        V2=VX**2+VY**2+VZ**2
        COLOUR=MIN(XRV**2/(ABS(XL2*V2)+SMALL2),1.0)
      ELSE
        COLOUR=MIN(-XNL/SQRT(ABS(XL2*XN2)+SMALL2),1.0)
      ENDIF
      END SUBROUTINE SB2SR4

      SUBROUTINE SBPLNT(EYE,NV,VERT,IC1,IC2,LIGHT,ITRANS)
C     ---------------------------------------------------
C
      REAL             EYE(*),VERT(3,*),LIGHT(*)
C
      REAL*8           XLNORM,ZZ,DZZ
      REAL             XW(400),YW(400)
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
!




C      This subroutine plots a diffusively-lit, semi-transparent, 
C    coloured plane; the use must ensure that all the verticies lie in a
C    flat plane, and that the bounding polygon be convex (so that the 
C    angle at any vertex <= 180 degs). All (x,y,z) values are taken to 
C    be given in world coordinates. The z-component of the eye-poisition
C    should be positive and that of the vertices should be negative; the
C    viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    NV       R*4    I       -      No. of verticies (>=3).
C    VERT     R*4    I     3 x NV   (x,y,z) coordinate of verticies.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    ITRANS   I*4    I       -      Level of transparency:
C                                        1 = 25%; 2 = 50%; 3 = 75%.
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBLIN1     Calculates the projection of (x,y,z) on viewing screen.
C
C History
C   D. S. Sivia      21 Aug 1995  Initial release.
C   D. S. Sivia      24 Oct 1997 "Safe-guarded" some rounding errors.
C-----------------------------------------------------------------------
C
C Carry out some initial checks and calculate the coordinates of the 
C projected polygon.
C
      SMALL=1.0E-10
      SMALL2=SMALL**2
      IF (EYE(3).LE.SMALL) RETURN
      IF (NV.LT.3 .OR. NV.GT.400) RETURN
      DO 10 I=1,NV
  10    IF (VERT(3,I).GE.0.0) RETURN
      XMIN=+1.0E20
      XMAX=-1.0E20
      YMIN=+1.0E20
      YMAX=-1.0E20
      DO 20 I=1,NV
        CALL SBLIN1(EYE,VERT(1,I),VERT(2,I),VERT(3,I),XW(I),YW(I))
        IF (XW(I).LT.XMIN) THEN
          XMIN=XW(I)
          ILEFT=I
        ENDIF
        IF (YW(I).LT.YMIN) THEN
          YMIN=YW(I)
          JBOTOM=I
        ENDIF
        XMAX=MAX(XW(I),XMAX)
        YMAX=MAX(YW(I),YMAX)
  20  CONTINUE
      IF (XMIN.GE.XTRC .OR. XMAX.LE.XBLC) RETURN
      IF (YMIN.GE.YTRC .OR. YMAX.LE.YBLC) RETURN
C
C Find the outward normal seen by the eye, and activate the appropriate 
C colour.
C
      AX=VERT(1,2)-VERT(1,1)
      AY=VERT(2,2)-VERT(2,1)
      AZ=VERT(3,2)-VERT(3,1)
      BX=VERT(1,1)-VERT(1,NV)
      BY=VERT(2,1)-VERT(2,NV)
      BZ=VERT(3,1)-VERT(3,NV)
      XN=BY*AZ-AY*BZ
      YN=BZ*AX-AZ*BX
      ZN=BX*AY-AX*BY
      TEN=XN*(EYE(1)-VERT(1,1))+YN*(EYE(2)-VERT(2,1))
     *   +ZN*(EYE(3)-VERT(3,1))
      COLOUR=FLOAT(IC1)
      NC=IC2-IC1
      IF (NC.GT.0) THEN
        TNL=XN*LIGHT(1)+YN*LIGHT(2)+ZN*LIGHT(3)
        IF (TEN.LT.0.0) TNL=-TNL
        COSDIF=0.0
        IF (TNL.LT.0.0) THEN
          TN2=XN**2+YN**2+ZN**2
          TL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
          COSDIF=MIN(-TNL/SQRT(TN2*TL2+SMALL2),1.0)
        ENDIF
        COLOUR=COLOUR+COSDIF*FLOAT(NC)
      ENDIF
C
C Plot the projected polygon.
C
      ITLEVL=MAX(MIN(ITRANS,3),1)
      XLNORM=DBLE(EYE(3))*DBLE(TEN)
      EYENRM=XN*EYE(1)+YN*EYE(2)+ZN*EYE(3)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      DYJ=1.0/DY
      DXI=1.0/DX
      SAFER=0.0001
      IF ((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
        JYMIN=INT((YMIN-YBLC)*DY)+2
        JYMAX=MIN(INT((YMAX-YBLC)*DY)+1,NYP)
        IF (JYMIN.GT.JYMAX) RETURN
        YJ=YBLC+(FLOAT(JYMIN-1)+SAFER)*DYJ
        NVL2=JBOTOM
        NVR2=JBOTOM
        J1=JYMIN
        DO 50 IVERT=1,NV
          IF (YJ.GT.YW(NVL2)) THEN
   1        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVL2)) GOTO 1
            YDIFL=YW(NVL2)-YW(NVL1)
            IF (ABS(YDIFL).LT.SMALL) YDIFL=SMALL
            GRADL=(XW(NVL2)-XW(NVL1))/YDIFL
          ENDIF
          IF (YJ.GT.YW(NVR2)) THEN
   2        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVR2)) GOTO 2
            YDIFR=YW(NVR2)-YW(NVR1)
            IF (ABS(YDIFR).LT.SMALL) YDIFR=SMALL
            GRADR=(XW(NVR2)-XW(NVR1))/YDIFR
          ENDIF
          IF (YW(NVL2).LT.YW(NVR2)) THEN
            J2=MIN(INT((YW(NVL2)-YBLC)*DY)+1,JYMAX)
          ELSE
            J2=MIN(INT((YW(NVR2)-YBLC)*DY)+1,JYMAX)
          ENDIF
          DO 40 J=J1,J2
            IF (J.GE.1) THEN
              JTEST=MOD(J,2)
              IF (ITLEVL.EQ.3 .AND. JTEST.EQ.1) GOTO 39
              XL=XW(NVL1)+GRADL*(YJ-YW(NVL1))
              XR=XW(NVR1)+GRADR*(YJ-YW(NVR1))
              ISTEP=1
              IX1=MAX(INT((XL-XBLC)*DX)+2,1)
              IX2=MIN(INT((XR-XBLC)*DX)+1,NXP)
              IF (IX1.GT.IX2) THEN
                ISTEP=-1
                IX1=MIN(IX1-1,NXP)
                IX2=MAX(IX2+1,1)
              ENDIF
              DZZ=DBLE(FLOAT(ISTEP)*DXI*XN)
              ZZ=DBLE(EYENRM-(XBLC+FLOAT(IX1-1)*DXI)*XN-YJ*YN)
              K=(J-1)*NXP+IX1
              DO 30 I=IX1,IX2,ISTEP
                ITEST=MOD(I,2)
                IF (ITLEVL.EQ.1) THEN
                  IF ((ITEST+JTEST).EQ.0) GOTO 29
                ELSEIF (ITLEVL.EQ.2) THEN
                  IF ((ITEST+JTEST).EQ.1) GOTO 29
                ELSE
                  IF (ITEST.EQ.1) GOTO 29
                ENDIF
                Z=EYE(3)-SNGL(XLNORM/ZZ)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  SBBUFF(KSTART+K)=COLOUR
                ENDIF
  29            ZZ=ZZ-DZZ
                K=K+ISTEP
  30          CONTINUE
            ENDIF
  39        YJ=YJ+DYJ
  40      CONTINUE
          J1=J2+1
          IF (J1.GT.JYMAX) RETURN
  50    CONTINUE
      ELSE
        IXMIN=INT((XMIN-XBLC)*DX)+2
        IXMAX=MIN(INT((XMAX-XBLC)*DX)+1,NXP)
        IF (IXMIN.GT.IXMAX) RETURN
        XI=XBLC+(FLOAT(IXMIN-1)+SAFER)*DXI
        NVL2=ILEFT
        NVR2=ILEFT
        I1=IXMIN
        DO 80 IVERT=1,NV
          IF (XI.GT.XW(NVL2)) THEN
   3        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVL2)) GOTO 3
            XDIFL=XW(NVL2)-XW(NVL1)
            IF (ABS(XDIFL).LT.SMALL) XDIFL=SMALL
            GRADL=(YW(NVL2)-YW(NVL1))/XDIFL
          ENDIF
          IF (XI.GT.XW(NVR2)) THEN
   4        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVR2)) GOTO 4
            XDIFR=XW(NVR2)-XW(NVR1)
            IF (ABS(XDIFR).LT.SMALL) XDIFR=SMALL
            GRADR=(YW(NVR2)-YW(NVR1))/XDIFR
          ENDIF
          IF (XW(NVL2).LT.XW(NVR2)) THEN
            I2=MIN(INT((XW(NVL2)-XBLC)*DX)+1,IXMAX)
          ELSE
            I2=MIN(INT((XW(NVR2)-XBLC)*DX)+1,IXMAX)
          ENDIF
          DO 70 I=I1,I2
            IF (I.GE.1) THEN
              ITEST=MOD(I,2)
              IF (ITLEVL.EQ.3 .AND. ITEST.EQ.1) GOTO 69
              YL=YW(NVL1)+GRADL*(XI-XW(NVL1))
              YR=YW(NVR1)+GRADR*(XI-XW(NVR1))
              ISTEP=1
              JY1=MAX(INT((YL-YBLC)*DY)+2,1)
              JY2=MIN(INT((YR-YBLC)*DY)+1,NYP)
              IF (JY1.GT.JY2) THEN
                ISTEP=-1
                JY1=MIN(JY1-1,NYP)
                JY2=MAX(JY2+1,1)
              ENDIF
              DZZ=DBLE(FLOAT(ISTEP)*DYJ*YN)
              ZZ=DBLE(EYENRM-(YBLC+FLOAT(JY1-1)*DYJ)*YN-XI*XN)
              K=(JY1-1)*NXP+I
              KSTEP=ISTEP*NXP
              DO 60 J=JY1,JY2,ISTEP
                JTEST=MOD(J,2)
                IF (ITLEVL.EQ.1) THEN
                  IF ((ITEST+JTEST).EQ.0) GOTO 59
                ELSEIF (ITLEVL.EQ.2) THEN
                  IF ((ITEST+JTEST).EQ.1) GOTO 59
                ELSE
                  IF (JTEST.EQ.1) GOTO 59
                ENDIF
                Z=EYE(3)-SNGL(XLNORM/ZZ)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  SBBUFF(KSTART+K)=COLOUR
                ENDIF
  59            ZZ=ZZ-DZZ
                K=K+KSTEP
  60          CONTINUE
            ENDIF
  69        XI=XI+DXI
  70      CONTINUE
          I1=I2+1
          IF (I1.GT.IXMAX) RETURN
  80    CONTINUE
      ENDIF
      END SUBROUTINE SBPLNT

      SUBROUTINE SBELIP(EYE,CENTRE,PAXES,IC1,IC2,LIGHT,LSHINE,ICLINE,
     *                  ANGLIN,X0,Y0,R0)
C     ---------------------------------------------------------------
C
      REAL            EYE(*),CENTRE(*),PAXES(3,*),LIGHT(*)
      LOGICAL         LSHINE
C
      REAL            SURF(3),EVEC(3,3),ENRM(3)
      REAL*8          DXE,DYE,DZE,DZE2,DWZE,DXH,DYH,DUXH,DVYH
      REAL*8          ALFA,BETA,GAMA,QXX,QYY,QZZ,QXY,QYZ,QZX,U,V,W,XMU
      REAL*8          DBL0,DBL1,DBL2,A,B,C,DET,Q,DX0H,DY0H,DSMALL
      REAL*8          XL0,XL1,HYP,SINPHI,COSPHI,R1,R2
      LOGICAL         LPS,LCOLOR
      COMMON /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
!




C      This subroutine plots a shiny or matt coloured elliptical ball. 
C    All (x,y,z) values are taken to be given in world coordinates. The 
C    z-component of the eye-poisition should be positive and that of
C    the ball-centre should be negative (< -radius); the viewing-screen
C    is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    CENTRE   R*4    I       3      (x,y,z) coordinate of ball-centre.
C    PAXES    R*4    I     3 x 3    Principal axes of the elliposid.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    LSHINE   L*1    I       -      Shiny ball if .TRUE., else diffuse.
C    ICLINE   I*4    I       -      If >=0, colour index for lines on
C                                   surface of ellipsoid.
C    ANGLIN   R*4    I       -      Width of lines: +/- degs.
C    X0,Y0    R*4    O       -      Centre of projected ball.
C    R0       R*4    O       -      Average radius of projected ball.
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBEGLS     Works out colour-shade for surface of ellipsoid.
C
C History
C   D. S. Sivia       8 Sep 1995  Initial release.
C   D. S. Sivia      29 Sep 1995  Pass down angular-width of lines.
C-----------------------------------------------------------------------
C
C Set up standard ellipsoid and carry out some initial checks.
C
      SMALL=1.0E-20
      DSMALL=DBLE(SMALL)
      IF (EYE(3).LE.0.0) RETURN
      EVAL1=PAXES(1,1)**2+PAXES(2,1)**2+PAXES(3,1)**2
      EVAL2=PAXES(1,2)**2+PAXES(2,2)**2+PAXES(3,2)**2
      EVAL3=PAXES(1,3)**2+PAXES(2,3)**2+PAXES(3,3)**2
      IF (EVAL1.LT.SMALL .OR. EVAL2.LT.SMALL .OR. EVAL3.LT.SMALL) RETURN
      EVAL1=1.0/EVAL1
      EVAL2=1.0/EVAL2
      EVAL3=1.0/EVAL3
      ENRM(1)=SQRT(EVAL1)
      ENRM(2)=SQRT(EVAL2)
      ENRM(3)=SQRT(EVAL3)
      DO 20 J=1,3
        DO 10 I=1,3
  10      EVEC(I,J)=PAXES(I,J)*ENRM(J)
  20  CONTINUE
      DOTSUM=EVEC(1,1)*EVEC(1,2)+EVEC(2,1)*EVEC(2,2)+EVEC(3,1)*EVEC(3,2)
     *      +EVEC(1,2)*EVEC(1,3)+EVEC(2,2)*EVEC(2,3)+EVEC(3,2)*EVEC(3,3)
     *      +EVEC(1,3)*EVEC(1,1)+EVEC(2,3)*EVEC(2,1)+EVEC(3,3)*EVEC(3,1)
      IF (DOTSUM.GT.0.001) RETURN
      QXX=DBLE(EVAL1*EVEC(1,1)**2+EVAL2*EVEC(1,2)**2+EVAL3*EVEC(1,3)**2)
      QYY=DBLE(EVAL1*EVEC(2,1)**2+EVAL2*EVEC(2,2)**2+EVAL3*EVEC(2,3)**2)
      QZZ=DBLE(EVAL1*EVEC(3,1)**2+EVAL2*EVEC(3,2)**2+EVAL3*EVEC(3,3)**2)
      QXY=DBLE(EVAL1*EVEC(1,1)*EVEC(2,1)+EVAL2*EVEC(1,2)*EVEC(2,2)
     *        +EVAL3*EVEC(1,3)*EVEC(2,3))
      QYZ=DBLE(EVAL1*EVEC(2,1)*EVEC(3,1)+EVAL2*EVEC(2,2)*EVEC(3,2)
     *        +EVAL3*EVEC(2,3)*EVEC(3,3))
      QZX=DBLE(EVAL1*EVEC(3,1)*EVEC(1,1)+EVAL2*EVEC(3,2)*EVEC(1,2)
     *        +EVAL3*EVEC(3,3)*EVEC(1,3))
      DBL1=(QYZ*QXX-QXY*QZX)/(QXY**2-QXX*QYY)
      DBL2=(QYZ*QXY-QYY*QZX)/(QXX*QYY-QXY**2)
      A=QXX*DBL2**2+QYY*DBL1**2+QZZ
      B=QXY*DBL1*DBL2+QYZ*DBL1+QZX*DBL2
      ZMAX=CENTRE(3)+SNGL((SQRT(B**2+A)-B)/A)
      IF (ZMAX.GT.-SMALL) RETURN
C
C Calculate some useful parameters.
C
      DXE=DBLE(EYE(1))
      DYE=DBLE(EYE(2))
      DZE=DBLE(EYE(3))
      DZE2=DZE**2
      ALFA=DBLE(EYE(1)-CENTRE(1))
      BETA=DBLE(EYE(2)-CENTRE(2))
      GAMA=DBLE(EYE(3)-CENTRE(3))
      U=ALFA*QXX+BETA*QXY+GAMA*QZX
      V=ALFA*QXY+BETA*QYY+GAMA*QYZ
      W=ALFA*QZX+BETA*QYZ+GAMA*QZZ
      XMU=QXX*ALFA**2+QYY*BETA**2+QZZ*GAMA**2+2.0D0*(ALFA*BETA*QXY
     *    +BETA*GAMA*QYZ+GAMA*ALFA*QZX)-1.0D0
      A=XMU*QXX-U**2
      B=XMU*QXY-U*V
      C=XMU*QYY-V**2
      DET=ABS(A*C-B**2)+DSMALL
      DBL1=DZE*(XMU*QZX-U*W)
      DBL2=DZE*(XMU*QYZ-V*W)
      DX0H=(DBL1*C-DBL2*B)/DET
      DY0H=(DBL2*A-DBL1*B)/DET
      Q=DZE2*(W**2-XMU*QZZ)+A*DX0H**2+2.0D0*B*DX0H*DY0H+C*DY0H**2
      X0H=SNGL(DX0H)
      Y0H=SNGL(DY0H)
      DX=(XTRC-XBLC)/FLOAT(NXP-1)
      DY=(YTRC-YBLC)/FLOAT(NYP-1)
      XDIF=SNGL(SQRT(ABS(C*Q/DET)+DSMALL))
      XMIN=X0H-XDIF+EYE(1)
      XMAX=X0H+XDIF+EYE(1)
      IXMIN=INT((XMIN-XBLC)/DX)+2
      IXMAX=INT((XMAX-XBLC)/DX)+1
      IF (IXMIN.GT.NXP .OR. IXMAX.LT.1) RETURN
      YDIF=(SQRT(ABS(A*Q/DET)+DSMALL))
      YMIN=Y0H-YDIF+EYE(2)
      YMAX=Y0H+YDIF+EYE(2)
      JYMIN=INT((YMIN-YBLC)/DY)+2
      JYMAX=INT((YMAX-YBLC)/DY)+1
      IF (JYMIN.GT.NYP .OR. JYMAX.LT.1) RETURN
      IF (JYMIN.LT.1) JYMIN=1
      IF (JYMAX.GT.NYP) JYMAX=NYP
      X0=X0H+EYE(1)
      Y0=Y0H+EYE(2)
      COREL=SNGL(SQRT(ABS((B*B)/(A*C))+DSMALL))
      IF (COREL.GT.0.0001) THEN
        XL0=(A+C)/2.0D0
        XL1=XL0-SQRT(ABS(XL0*XL0-DET)+DSMALL)
        HYP=SQRT((XL1-A)**2+B**2+DSMALL)
        SINPHI=(XL1-A)/HYP
        COSPHI=B/HYP
      ELSE
        SINPHI=0.0D0
        COSPHI=1.0D0
      ENDIF
      R1=SQRT(Q/(A*COSPHI*COSPHI+SINPHI*(C*SINPHI+2.0*B*COSPHI)))
      R2=SQRT(Q/(A*SINPHI*SINPHI+COSPHI*(C*COSPHI-2.0*B*SINPHI)))
      R0=SNGL((R1+R2)/2.0D0)
C
C Fill the inside of the projected ellipse with the right colours.
C
      COSMIN=MAX(SIN(ABS(ANGLIN/57.29578)),0.0001)
      NC=IC2-IC1
      COL0=FLOAT(IC1)
      COLSCL=FLOAT(NC)
      COLINE=FLOAT(ICLINE)
      XL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
      YH=YBLC+DY*FLOAT(JYMIN-1)-EYE(2)
      DWZE=W*DZE
      BA=SNGL(B/A)
      DO 40 JY=JYMIN,JYMAX
        DYH=DBLE(YH)
        DVYH=V*DYH
        DBL0=QZZ*DZE2+DYH*(QYY*DYH-2.0D0*QYZ*DZE)
        YDIF=YH-Y0H
        XDIF=SNGL(SQRT(ABS(A*Q-DET*DBLE(YDIF**2))+DSMALL)/A)
        XMIN=X0H-BA*YDIF-XDIF+EYE(1)
        XMAX=X0H-BA*YDIF+XDIF+EYE(1)
        IXMIN=INT((XMIN-XBLC)/DX)+2
        IXMAX=INT((XMAX-XBLC)/DX)+1
        IF (IXMIN.LE.NXP .AND. IXMAX.GE.1) THEN
          IF (IXMIN.LT.1) IXMIN=1
          IF (IXMAX.GT.NXP) IXMAX=NXP
          XH=XBLC+DX*FLOAT(IXMIN-1)-EYE(1)
          K=(JY-1)*NXP+IXMIN
          DO 30 IX=IXMIN,IXMAX
            IF (ZMAX.GT.SBBUFF(K)) THEN
              DXH=DBLE(XH)
              DUXH=U*DXH
              DBL1=DUXH+DVYH-DWZE
              DBL2=DBL0+DXH*(QXX*DXH+2.0D0*(QXY*DYH-QZX*DZE))
              XLM=SNGL((-DBL1-SQRT(ABS(DBL1**2-XMU*DBL2)+DSMALL))/DBL2)
              SURF(3)=EYE(3)*(1.0-XLM)
              IF (SURF(3).GT.SBBUFF(K)) THEN
                SBBUFF(K)=SURF(3)
                IF (NC.EQ.0) THEN
                  SBBUFF(KSTART+K)=COL0
                ELSE
                  SURF(2)=EYE(2)+YH*XLM
                  SURF(1)=EYE(1)+XH*XLM
                  CALL SBEGLS(EYE,CENTRE,LIGHT,SURF,XL2,QXX,QYY,QZZ,QXY,
     *              QYZ,QZX,SMALL,LSHINE,COLOUR,EVEC,ICLINE,COSMIN)
                  IF (COLOUR.GT.2.0) THEN
                    SBBUFF(KSTART+K)=COLINE
                  ELSE
                    SBBUFF(KSTART+K)=COL0+COLOUR*COLSCL
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            K=K+1
            XH=XH+DX
  30      CONTINUE
        ENDIF
        YH=YH+DY
  40  CONTINUE
      END SUBROUTINE SBELIP

      SUBROUTINE SBEGLS(EYE,CENTRE,LIGHT,SURF,XL2,QXX,QYY,QZZ,QXY,QYZ,
     *                  QZX,SMALL,LSHINE,COLOUR,EVEC,ICLINE,COSANG)
C     ----------------------------------------------------------------
C
C Support subroutine for SBELIP, to work out colour-shade.
C
      REAL    EYE(*),CENTRE(*),LIGHT(*),SURF(*),EVEC(3,*)
      REAL*8  QXX,QYY,QZZ,QXY,QYZ,QZX,DDX,DDY,DDZ
      LOGICAL LSHINE
      REAL    NORMAL(3),REFLEC(3),VIEW(3)
C
      COLOUR=0.0
      DX=SURF(1)-CENTRE(1)
      DY=SURF(2)-CENTRE(2)
      DZ=SURF(3)-CENTRE(3)
      IF (ICLINE.GE.0) THEN
        SNORM=1.0/SQRT(DX**2+DY**2+DZ**2)
        COS1=ABS(SNORM*(DX*EVEC(1,1)+DY*EVEC(2,1)+DZ*EVEC(3,1)))
        COS2=ABS(SNORM*(DX*EVEC(1,2)+DY*EVEC(2,2)+DZ*EVEC(3,2)))
        COS3=ABS(SNORM*(DX*EVEC(1,3)+DY*EVEC(2,3)+DZ*EVEC(3,3)))
        COSMIN=MIN(COS1,COS2,COS3)
        IF (COSMIN.LT.COSANG) THEN
          COLOUR=10.0
          RETURN
        ENDIF
      ENDIF
      DDX=DBLE(DX)
      DDY=DBLE(DY)
      DDZ=DBLE(DZ)
      NORMAL(1)=SNGL(QXX*DDX+QXY*DDY+QZX*DDZ)
      NORMAL(2)=SNGL(QXY*DDX+QYY*DDY+QYZ*DDZ)
      NORMAL(3)=SNGL(QZX*DDX+QYZ*DDY+QZZ*DDZ)
      XN2=0.0
      XNL=0.0
      DO 10 I=1,3
        XNL=XNL+NORMAL(I)*LIGHT(I)
        XN2=XN2+NORMAL(I)**2
  10  CONTINUE
      IF (XNL.GE.0.0) RETURN
      IF (LSHINE) THEN
        RFNORM=(XNL+XNL)/(XN2+SMALL)
        XRV=0.0
        DO 20 I=1,3
          VIEW(I)=EYE(I)-SURF(I)
          REFLEC(I)=LIGHT(I)-RFNORM*NORMAL(I)
          XRV=XRV+REFLEC(I)*VIEW(I)
  20    CONTINUE
        IF (XRV.LT.0.0) RETURN
        REF2=0.0
        VIEW2=0.0
        DO 30 I=1,3
          REF2=REF2+REFLEC(I)**2
          VIEW2=VIEW2+VIEW(I)**2
  30    CONTINUE
        COLOUR=MIN(XRV**2/(ABS(REF2*VIEW2)+SMALL),1.0)
      ELSE
        COLOUR=MIN(-XNL/SQRT(XN2*XL2+SMALL),1.0)
      ENDIF
      END SUBROUTINE SBEGLS

      SUBROUTINE SBTEXT(EYE,TEXT,ICOL,PIVOT,FJUST,ORIENT,SIZE)
C     --------------------------------------------------------
C
      use grpckg1inc
      REAL          EYE(*),PIVOT(*),ORIENT(3,2)
      CHARACTER*(*) TEXT
C
      REAL          END1(3),END2(3)
      INTEGER       SYMBOL(256),XYGRID(300),XYBASE,XYLEFT
      LOGICAL       LSTART,UNUSED
      SAVE          LSTART
      DATA          LSTART,SMALL,IJFLAG /.FALSE.,1.0E-10,-64/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Write a text string in 3-d perspective. All (x,y,z) values are 
C    taken to be given in world coordinates. The z-component of the 
C    eye-poisition should be positive and that of the text string should 
C    be negative; the viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    TEXT     C*1    I       *      The text string to be written.
C    ICOL     I*4    I       -      Colour index for text.
C    PIVOT    R*4    I       3      (x,y,z) coordinate of pivot point.
C    FJUST    R*4    I       -      Position of pivot along the text: 
C                                   0.0=left, 0.5=centre, 1.0=right.
C    ORIENT   R*4    I     3 x 2    (x,y,z) for X-length and Y-height
C                                   directions of the text.
C    SIZE     R*4    I       -      Height of the reference symbol "A".
C
C Globals
C    GRPCKG1.INC
C
C External Calls
!




C   SUBROUTINE   DESCRIPTION
C     GRSY00     Initialises font description.
C     GRSYDS     Decodes string into a list of Hershey symbol numbers.
C     GRSYXD     Obtains the polyline representation of a given symbol.
C     GRLEN      Calculates the length of the string.
C     SBLINE     Draws a 3-d line in perspective.
C
C History
C   D. S. Sivia      14 Sep 1995  Initial release.
C   D. S. Sivia      11 Oct 1995  Modified a J DO LOOP for a Pentium!
C-----------------------------------------------------------------------
C
C Carry out some initial checks.
C
      IF (EYE(3).LE.SMALL) RETURN
      IF (FJUST.LT.0.0 .OR. FJUST.GT.1.0) RETURN
      IF (SIZE.LT.SMALL) RETURN
      XLEN=SQRT(ORIENT(1,1)**2+ORIENT(2,1)**2+ORIENT(3,1)**2)
      YLEN=SQRT(ORIENT(1,2)**2+ORIENT(2,2)**2+ORIENT(3,2)**2)
      IF (XLEN.LT.SMALL .OR. YLEN.LT.SMALL) RETURN
      COSANG=(ORIENT(1,1)*ORIENT(1,2)+ORIENT(2,1)*ORIENT(2,2)+
     *        ORIENT(3,1)*ORIENT(3,2))/(XLEN*YLEN)
      IF (ABS(COSANG).GT.0.001) RETURN
      IF (.NOT.LSTART) THEN
        CALL GRSY00
        LSTART=.TRUE.
      ENDIF
      NCHAR=MIN(LEN(TEXT),256)
      DO 10 I=NCHAR,1,-1
  10    IF (TEXT(I:I).NE.' ') GOTO 1
   1  NCHMAX=I
      IF (NCHMAX.LT.1) RETURN
      DO 20 I=1,NCHMAX
  20    IF (TEXT(I:I).NE.' ') GOTO 2
   2  NCHMIN=I
C
C Calculate the parameters for the Hershey --> world coordinates 
C transformation.
C
      CALL GRLEN(TEXT(NCHMIN:NCHMAX),D)
      D=D*2.5/GRCFAC(GRCIDE)
      CALL GRSYDS(ISYMBA,NSYMBS,'A',1)
      CALL GRSYXD(ISYMBA,XYGRID,UNUSED)
      TSCL=SIZE/FLOAT(XYGRID(3)-XYGRID(2))
      XNORM=TSCL/XLEN
      YNORM=TSCL/YLEN
      XX=ORIENT(1,1)*XNORM
      XY=ORIENT(2,1)*XNORM
      XZ=ORIENT(3,1)*XNORM
      YX=ORIENT(1,2)*YNORM
      YY=ORIENT(2,2)*YNORM
      YZ=ORIENT(3,2)*YNORM
      X0=PIVOT(1)-FJUST*D*XX
      Y0=PIVOT(2)-FJUST*D*XY
      Z0=PIVOT(3)-FJUST*D*XZ
      Z1=PIVOT(3)+(1.0-FJUST)*D*XZ
      IF (Z0.GT.-SMALL .OR. Z1.GT.-SMALL) RETURN
C
C Write out the text string, character by character.
C
      DX=0.0
      DY=0.0
      DZ=0.0
      FNTBAS=0.0
      FNTFAC=1.0
      IFNTLV=0
      CALL GRSYDS(SYMBOL,NSYMBS,TEXT(NCHMIN:NCHMAX),1)
      DO 40 K=1,NSYMBS
        KSYMB=SYMBOL(K)
        IF (KSYMB.LT.0) THEN
          IF (KSYMB.EQ.-1) THEN
            IFNTLV=IFNTLV+1
            FNTBAS=FNTBAS+16.0*FNTFAC
            FNTFAC=0.75**ABS(IFNTLV)
          ELSEIF (KSYMB.EQ.-2) THEN
            IFNTLV=IFNTLV-1
            FNTFAC=0.75**ABS(IFNTLV)
            FNTBAS=FNTBAS-16.0*FNTFAC
          ELSEIF (KSYMB.EQ.-3) THEN
            X0=X0-DX
            Y0=Y0-DY
            Z0=Z0-DZ
          ENDIF
          GOTO 40
        ENDIF
        CALL GRSYXD(KSYMB,XYGRID,UNUSED)
        IF (.NOT. UNUSED) THEN
          XYBASE=XYGRID(2)
          XYLEFT=XYGRID(4)
          RLX=FLOAT(XYGRID(6)-XYLEFT)*FNTFAC
          RLY=FLOAT(XYGRID(7)-XYBASE)*FNTFAC+FNTBAS
          END1(1)=X0+RLX*XX+RLY*YX
          END1(2)=Y0+RLX*XY+RLY*YY
          END1(3)=Z0+RLX*XZ+RLY*YZ
          J=8
          DO 30 JJ=8,298,2
            IX=XYGRID(J)
            JY=XYGRID(J+1)
            IF (JY.EQ.IJFLAG) GOTO 3
            IF (IX.EQ.IJFLAG) THEN
              J=J+2
              RLX=FLOAT(XYGRID(J)-XYLEFT)*FNTFAC
              RLY=FLOAT(XYGRID(J+1)-XYBASE)*FNTFAC+FNTBAS
              END1(1)=X0+RLX*XX+RLY*YX
              END1(2)=Y0+RLX*XY+RLY*YY
              END1(3)=Z0+RLX*XZ+RLY*YZ
            ELSE
              RLX=FLOAT(IX-XYLEFT)*FNTFAC
              RLY=FLOAT(JY-XYBASE)*FNTFAC+FNTBAS
              END2(1)=X0+RLX*XX+RLY*YX
              END2(2)=Y0+RLX*XY+RLY*YY
              END2(3)=Z0+RLX*XZ+RLY*YZ
              CALL SBLINE(EYE,END1,END2,ICOL,.FALSE.)
              CALL SBRCOP(END2,END1,3)
            ENDIF
            J=J+2            
  30      CONTINUE
        ENDIF
  3     XL=FNTFAC*FLOAT(XYGRID(5)-XYLEFT)
        DX=XL*XX
        DY=XL*XY
        DZ=XL*XZ
        X0=X0+DX
        Y0=Y0+DY
        Z0=Z0+DZ
  40  CONTINUE
      END SUBROUTINE SBTEXT

      SUBROUTINE SBCPLN(EYE,LATICE,IC1,IC2,LIGHT,SLNORM,APOINT,ICEDGE,
     *                  ITRANS)
C     ----------------------------------------------------------------
C
      REAL    EYE(*),LATICE(3,*),LIGHT(*)
      REAL    SLNORM(*),APOINT(*)
C
      REAL    END1(3),END2(3),VERT1(3,12),VERT2(3,12)
      LOGICAL LVERT(12)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
!




C      This subroutine plots a diffusively-lit, semi-transparent, 
C    coloured plane through a unit cell. All (x,y,z) values are taken to
C    be given in world coordinates. The z-component of the eye-poisition 
C    should be positive and that of all the lattice-vertices should be 
C    negative; the viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    LATICE   R*4    I     3 x 4    (x,y,z) coordinates of the origin
C                                   and the a, b & C lattice-vertices.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    SLNORM   R*4    I       3      (x,y,z) direction of normal to plane.
C    APONIT   R*4    I       3      (x,y,z) coordinate of a point within
C                                   the plane.
C    ICEDGE   I*4    I       -      If >=0, it's the colour-index for
C                                   the boundary of the plane.
C    ITRANS   I*4    I       -      Level of transparency:
C                                     0 = 0%; 1 = 25%; 2 = 50%; 3 = 75%.
C
C Globals 
C    None
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBSLC1     Checks whether a side of the unit cell is cut by the
C                plane & calculates the coordinates of the intersection.
C     SBCPL1     Order the verticies of the polygon-of-intersection.
C     SBLIN1     Calculates the projection of (x,y,z) on viewing screen.
C     SBPLAN     Plots a coloured plane.
C     SBPLNT     Plots a semi-transparent coloured plane.
C
C History
C   D. S. Sivia      26 Sep 1995  Initial release.
C-----------------------------------------------------------------------
C
C Carry out some initial checks.
C
      SMALL=1.0E-10
      SMALL2=SMALL**2
      IF (EYE(3).LE.SMALL) RETURN
      SNRM=1.0/SQRT(SLNORM(1)**2+SLNORM(2)**2+SLNORM(3)**2+SMALL2)
      XN=SLNORM(1)*SNRM
      YN=SLNORM(2)*SNRM
      ZN=SLNORM(3)*SNRM
      XAE=EYE(1)-APOINT(1)
      YAE=EYE(2)-APOINT(2)
      ZAE=EYE(3)-APOINT(3)
      XLNORM=DBLE(XN*XAE)+DBLE(YN*YAE)+DBLE(ZN*ZAE)
      COSNRM=SNGL(XLNORM/SQRT(DBLE(XAE**2)+DBLE(YAE**2)+DBLE(ZAE**2)
     *      +DBLE(SMALL2)))
      IF (ABS(COSNRM).LT.0.001) RETURN
      DO 10 J=1,3
        DX=LATICE(1,J+1)-LATICE(1,1)
        DY=LATICE(2,J+1)-LATICE(2,1)
        DZ=LATICE(3,J+1)-LATICE(3,1)
        D2=DX**2+DY**2+DZ**2
        IF (D2.LT.SMALL2) RETURN
  10  CONTINUE
C
C Calculate the coordinates of the polygon-of-intersection between the
C the plane and the edges of the unit cell.
C
      NVERT=0
      II=0
      DO 20 L=1,12,4
        I=II+2
        J=MOD(II+1,3)+2
        K=MOD(II+2,3)+2
        CALL SBSLC1(LATICE,LATICE(1,K),XN,YN,ZN,APOINT,LVERT(L),
     *              VERT1(1,L),NVERT)
        END1(1)=LATICE(1,I)+LATICE(1,J)-LATICE(1,1)
        END1(2)=LATICE(2,I)+LATICE(2,J)-LATICE(2,1)
        END1(3)=LATICE(3,I)+LATICE(3,J)-LATICE(3,1)
        CALL SBSLC1(LATICE(1,I),END1,XN,YN,ZN,APOINT,LVERT(L+1),
     *              VERT1(1,L+1),NVERT)
        CALL SBSLC1(LATICE(1,J),END1,XN,YN,ZN,APOINT,LVERT(L+2),
     *              VERT1(1,L+2),NVERT)
        END2(1)=END1(1)+LATICE(1,K)-LATICE(1,1)
        END2(2)=END1(2)+LATICE(2,K)-LATICE(2,1)
        END2(3)=END1(3)+LATICE(3,K)-LATICE(3,1)
        CALL SBSLC1(END1,END2,XN,YN,ZN,APOINT,LVERT(L+3),
     *              VERT1(1,L+3),NVERT)
        II=II+1
  20  CONTINUE
      IF (NVERT.LT.3) RETURN
      CALL SBCPL1(EYE,LVERT,VERT1,VERT2,NVERT,XN,YN,ZN,ICEDGE,ZDLINE)
C
C Plot the plane.
C
      ITLEVL=MAX(MIN(ITRANS,3),0)
      IF (ITLEVL.EQ.0) THEN
        CALL SBPLAN(EYE,NVERT,VERT2,IC1,IC2,LIGHT)
      ELSE
        CALL SBPLNT(EYE,NVERT,VERT2,IC1,IC2,LIGHT,ITLEVL)
      ENDIF
      END SUBROUTINE SBCPLN

      SUBROUTINE SBCPL1(EYE,LVERT,VERT1,VERT2,NVERT,XN,YN,ZN,ICOL,ZDIF)
C     -----------------------------------------------------------------
C
      REAL    EYE(*),VERT1(3,*),VERT2(3,*),VERT3(3),ANGLE(12)
      INTEGER ISORT(12)
      LOGICAL LVERT(*)
C
      IV1=0
      XBAR=0.0
      YBAR=0.0
      ZBAR=0.0
      ZMIN=+1.0E20
      ZMAX=-1.0E20
      DO 10 K=1,12
        ZMIN=MIN(ZMIN,VERT1(3,K))
        ZMAX=MAX(ZMAX,VERT1(3,K))
        IF (LVERT(K)) THEN
          IF (IV1.LE.0) IV1=K
          XBAR=XBAR+VERT1(1,K)
          YBAR=YBAR+VERT1(2,K)
          ZBAR=ZBAR+VERT1(3,K)
        ENDIF
  10  CONTINUE
      ZDIF=(ZMAX-ZMIN)/5000.0
      XBAR=XBAR/FLOAT(NVERT)
      YBAR=YBAR/FLOAT(NVERT)
      ZBAR=ZBAR/FLOAT(NVERT)
      XREF=VERT1(1,IV1)-XBAR
      YREF=VERT1(2,IV1)-YBAR
      ZREF=VERT1(3,IV1)-ZBAR
      REFNRM=1.0/SQRT(XREF**2+YREF**2+ZREF**2+1.0E-20)
      XREF=XREF*REFNRM
      YREF=YREF*REFNRM
      ZREF=ZREF*REFNRM
      XNRM=YREF*ZN-YN*ZREF
      YNRM=ZREF*XN-ZN*XREF
      ZNRM=XREF*YN-XN*YREF
      J=1
      ANGLE(J)=0.0
      ISORT(J)=IV1
      CALL SBRCOP(VERT1(1,IV1),VERT2(1,J),3)
      DO 40 K=IV1+1,12
        IF (LVERT(K)) THEN
          J=J+1
          XVEC=VERT1(1,K)-XBAR
          YVEC=VERT1(2,K)-YBAR
          ZVEC=VERT1(3,K)-ZBAR
          X=XVEC*XREF+YVEC*YREF+ZVEC*ZREF
          Y=XVEC*XNRM+YVEC*YNRM+ZVEC*ZNRM
          ANGJ=ATAN2(Y,X)
          CALL SBRCOP(VERT1(1,K),VERT3,3)
          DO 20 I=1,J-1
  20        IF (ANGJ.LT.ANGLE(I)) GOTO 1
   1      II=I
           DO 30 I=J,II+1,-1
            CALL SBRCOP(VERT2(1,I-1),VERT2(1,I),3)
            ANGLE(I)=ANGLE(I-1)
            ISORT(I)=ISORT(I-1)
  30      CONTINUE
          CALL SBRCOP(VERT3,VERT2(1,II),3)
          ANGLE(II)=ANGJ
          ISORT(II)=K
        ENDIF
  40  CONTINUE
      IF (ICOL.GE.0.0) THEN
        DO 50 I=1,NVERT-1
          J=ISORT(I)
          K=ISORT(I+1)
          CALL SBLINE(EYE,VERT1(1,J),VERT1(1,K),ICOL,.FALSE.)
  50    CONTINUE
        CALL SBLINE(EYE,VERT1(1,K),VERT1(1,ISORT(1)),ICOL,.FALSE.)
      ENDIF
      END SUBROUTINE SBCPL1

      SUBROUTINE SBTBAL(EYE,CENTRE,RADIUS,IC1,IC2,LIGHT,LSHINE,X0,Y0,R0,
     *                  ITRANS)
C     ------------------------------------------------------------------
C
      REAL            EYE(*),CENTRE(*),LIGHT(*)
      LOGICAL         LSHINE

      REAL            SURF(3)
      REAL*8          ALFA,BETA,GAMA,XMU,A,B,C,DET,Q,DX0H,DY0H
      REAL*8          DZE,DZE2,DGAMZE,DBL1,DBL2,DSMALL
      REAL*8          XL0,XL1,HYP,SINPHI,COSPHI,R1,R2
      LOGICAL         LPS,LCOLOR
      COMMON /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots a semi-transparent shiny or matt coloured 
C    ball. All (x,y,z) values are taken to be given in world coordinates.
C    The z-component of the eye-poisition should be positive and that of
C    the ball-centre should be negative (< -radius); the viewing-screen
C    is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    CENTRE   R*4    I       3      (x,y,z) coordinate of ball-centre.
C    RADIUS   R*4    I       -      Radius of ball.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    LSHINE   L*1    I       -      Shiny ball if .TRUE., else diffuse.
C    X0,Y0    R*4    O       -      Centre of projected ball.
C    R0       R*4    O       -      Average radius of projected ball.
C    ITRANS   I*4    I       -      Level of transparency:
C                                        1 = 25%; 2 = 50%; 3 = 75%.
C
C Globals 
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBGLOS     Works out colour-shade for surface of ball.
C
C History
C   D. S. Sivia       8 Jul 1996  Initial release.
C-----------------------------------------------------------------------
C
C Some initial checks.
C
      SMALL=1.0E-20
      DSMALL=DBLE(SMALL)
      IF (EYE(3).LE.0.0) RETURN
      IF (RADIUS.LE.0.0) RETURN
      IF (CENTRE(3).GT.-RADIUS) RETURN
C
C Calculate parameters of projected ellipse.
C
      DZE=DBLE(EYE(3))
      DZE2=DZE**2
      ALFA=DBLE(EYE(1)-CENTRE(1))
      BETA=DBLE(EYE(2)-CENTRE(2))
      GAMA=DBLE(EYE(3)-CENTRE(3))
      XMU=DBLE(RADIUS**2)-(ALFA**2+BETA**2+GAMA**2)
      A=XMU*(XMU+ALFA**2)
      B=XMU*ALFA*BETA
      C=XMU*(XMU+BETA**2)
      DET=ABS(A*C-B**2)+DSMALL
      DX0H=GAMA*XMU*DZE*(ALFA*C-BETA*B)/DET
      DY0H=GAMA*XMU*DZE*(BETA*A-ALFA*B)/DET
      Q=A*DX0H**2+2.0*B*DX0H*DY0H+C*DY0H**2-XMU*(XMU+GAMA**2)*DZE2
      X0H=SNGL(DX0H)
      Y0H=SNGL(DY0H)
      DX=(XTRC-XBLC)/FLOAT(NXP-1)
      DY=(YTRC-YBLC)/FLOAT(NYP-1)
      XDIF=SNGL(SQRT(ABS(C*Q/DET)+DSMALL))
      XMIN=X0H-XDIF+EYE(1)
      XMAX=X0H+XDIF+EYE(1)
      IXMIN=INT((XMIN-XBLC)/DX)+2
      IXMAX=INT((XMAX-XBLC)/DX)+1
      IF (IXMIN.GT.NXP .OR. IXMAX.LT.1) RETURN
      YDIF=(SQRT(ABS(A*Q/DET)+DSMALL))
      YMIN=Y0H-YDIF+EYE(2)
      YMAX=Y0H+YDIF+EYE(2)
      JYMIN=INT((YMIN-YBLC)/DY)+2
      JYMAX=INT((YMAX-YBLC)/DY)+1
      IF (JYMIN.GT.NYP .OR. JYMAX.LT.1) RETURN
      IF (JYMIN.LT.1) JYMIN=1
      IF (JYMAX.GT.NYP) JYMAX=NYP
      ZMAX=CENTRE(3)+RADIUS
      X0=X0H+EYE(1)
      Y0=Y0H+EYE(2)
      COREL=SNGL(SQRT(ABS((B*B)/(A*C))+DSMALL))
      IF (COREL.GT.0.0001) THEN
        XL0=(A+C)/2.0D0
        XL1=XL0-SQRT(ABS(XL0*XL0-DET)+DSMALL)
        HYP=SQRT((XL1-A)**2+B**2+DSMALL)
        SINPHI=(XL1-A)/HYP
        COSPHI=B/HYP
      ELSE
        SINPHI=0.0D0
        COSPHI=1.0D0
      ENDIF
      R1=SQRT(Q/(A*COSPHI*COSPHI+SINPHI*(C*SINPHI+2.0*B*COSPHI)))
      R2=SQRT(Q/(A*SINPHI*SINPHI+COSPHI*(C*COSPHI-2.0*B*SINPHI)))
      R0=SNGL((R1+R2)/2.0D0)
C
C Fill the inside of the projected ellipse with the right colours.
C
      ITLEVL=MAX(MIN(ITRANS,3),1)
      NC=IC2-IC1
      COL0=FLOAT(IC1)
      COLSCL=FLOAT(NC)
      XL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
      XN2=RADIUS**2
      XNL2=1.0/SQRT(XN2*XL2+SMALL)
      XN3=1.0/(XN2+SMALL)
      YH=YBLC+DY*FLOAT(JYMIN-1)-EYE(2)
      DGAMZE=GAMA*DZE
      BA=SNGL(B/A)
      DO 20 JY=JYMIN,JYMAX
        JTEST=MOD(JY,2)
        IF (ITLEVL.EQ.3 .AND. JTEST.EQ.1) GOTO 19
        YH2=YH**2
        BETAYH=BETA*YH
        YDIF=YH-Y0H
        XDIF=SNGL(SQRT(ABS(A*Q-DET*DBLE(YDIF**2))+DSMALL)/A)
        XMIN=X0H-BA*YDIF-XDIF+EYE(1)
        XMAX=X0H-BA*YDIF+XDIF+EYE(1)
        IXMIN=INT((XMIN-XBLC)/DX)+2
        IXMAX=INT((XMAX-XBLC)/DX)+1
        IF (IXMIN.LE.NXP .AND. IXMAX.GE.1) THEN
          IF (IXMIN.LT.1) IXMIN=1
          IF (IXMAX.GT.NXP) IXMAX=NXP
          XH=XBLC+DX*FLOAT(IXMIN-1)-EYE(1)
          K=(JY-1)*NXP+IXMIN
          DO 10 IX=IXMIN,IXMAX
            ITEST=MOD(IX,2)
            IF (ITLEVL.EQ.1) THEN
              IF ((ITEST+JTEST).EQ.0) GOTO 9
            ELSEIF (ITLEVL.EQ.2) THEN
              IF ((ITEST+JTEST).EQ.1) GOTO 9
            ELSE
              IF (ITEST.EQ.1) GOTO 9
            ENDIF
            IF (ZMAX.GT.SBBUFF(K)) THEN
              XH2=XH**2
              ALFAXH=ALFA*XH
              DBL1=DBLE(ALFAXH+BETAYH)-DGAMZE
              DBL2=DBLE(XH2+YH2)+DZE2
              XLM=SNGL((-DBL1-SQRT(ABS(DBL1**2+XMU*DBL2)+DSMALL))/DBL2)
              SURF(3)=EYE(3)*(1.0-XLM)
              IF (SURF(3).GT.SBBUFF(K)) THEN
                SBBUFF(K)=SURF(3)
                IF (NC.EQ.0) THEN
                  SBBUFF(KSTART+K)=COL0
                ELSE
                  SURF(2)=EYE(2)+YH*XLM
                  SURF(1)=EYE(1)+XH*XLM
                  CALL SBGLOS(EYE,CENTRE,LIGHT,SURF,XNL2,XN3,SMALL,
     *                        LSHINE,COLOUR)
                  SBBUFF(KSTART+K)=COL0+COLOUR*COLSCL
                ENDIF
              ENDIF
            ENDIF
   9        K=K+1
            XH=XH+DX
  10      CONTINUE
        ENDIF
  19    YH=YH+DY
  20  CONTINUE
      END SUBROUTINE SBTBAL

      SUBROUTINE SBTSUR(EYE,LATICE,DENS,N1,N2,N3,DSURF,IC1,IC2,LIGHT,
     *                  LSHINE,ITRANS)
C     ---------------------------------------------------------------
C
      REAL             EYE(*),LATICE(3,*),DENS(0:N1,0:N2,0:N3),LIGHT(*)
      LOGICAL          LSHINE
C
      REAL*8           XLNORM,PNEYE(3),DSMAL2
      REAL             BAS(3,3),MTRX
      REAL             XYZ(3),DXYZ(3,3),FRCXYZ(12),DDXYZ(3,12,2)
      REAL             DLOCAL(8),VERT(3,12),GRDSCL(3)
      LOGICAL          LPS,LCOLOR,LEMPTY
      INTEGER          IVERT(8)
      COMMON  /SRFCOM/ GRDCUB(3,8),MTRX(3,3),ORIG(3),XL2,COL0,COLSCL
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
!




C      This subroutine plots a semi-transparent iso-surface through a 
C    unit-cell of density. All (x,y,z) values are taken to be given in 
C    world coordinates. The z-component of the eye-poisition should be 
C    positive and that of all the lattice-vertices should be negative; 
C    the viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    LATICE   R*4    I     3 x 4    (x,y,z) coordinates of the origin
C                                   and the a, b & C lattice-vertices.
C    DENS     R*4    I     (N1+1)   The density at regular points within
C                        x (N2+1)   the unit cell, wrapped around so
C                        x (N3+1)   that DENS(0,J,K)=DENS(N1,J,K) etc..
C    N1,N2,N3 I*4    I       -      The dimensions of the unit-cell grid.
C    DSURF    R*4    I       -      Density for the iso-surface.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    LSHINE   L*1    I       -      Shiny surface if TRUE, else diffuse.
C    ITRANS   I*4    I       -      Level of transparency:
C                                        1 = 25%; 2 = 50%; 3 = 75%.
C
C Globals 
C    SFTBUF
C    SRFCOM
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBSRF0     A quick check in case there are no iso-surafces.
C     SBSRF1     Anlayses a 2-d box in the surface of the unit-cell.
C     SBTSF2     Paints a 2-d box in the surface of the unit-cell.
C     SBSRF3     Analyses a 3-d box within the unit-cell.
C     SBSRF4     Initialises the gradients for a 3-d box.
C     SBTSF5     Breaks up the iso-surface in a 3-d box into triangles.
C     SBTSF6     Paints a triangular patch of a semi-transp iso-surface.
C
C History
C   D. S. Sivia       9 Jul 1996  Initial release.
C   D. S. Sivia      24 Oct 1997 "Safe-guarded" some rounding errors.
C-----------------------------------------------------------------------
C
C Carry out some initial checks.
C
      SMALL=1.0E-10
      SMALL2=SMALL**2
      DSMAL2=DBLE(SMALL2)
      IF (EYE(3).LE.SMALL) RETURN
      IF (N1.LT.1 .OR. N2.LT.1 .OR. N3.LT.1) RETURN
      XL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
      IF (XL2.LT.SMALL) RETURN
      IF (LATICE(3,1).GE.0.0) RETURN
      ZFAR=LATICE(3,2)+LATICE(3,3)+LATICE(3,4)-2.0*LATICE(3,1)
      IF (ZFAR.GE.0.0) RETURN
      DO 10 J=1,3
        IF (LATICE(3,J+1).GE.0.0) RETURN
        BAS(1,J)=LATICE(1,J+1)-LATICE(1,1)
        BAS(2,J)=LATICE(2,J+1)-LATICE(2,1)
        BAS(3,J)=LATICE(3,J+1)-LATICE(3,1)
        IF ((ZFAR-BAS(3,J)).GE.0.0) RETURN
        BAS2J=BAS(1,J)**2+BAS(2,J)**2+BAS(3,J)**2
        IF (BAS2J.LT.SMALL2) RETURN
  10  CONTINUE
      NTOT=(N1+1)*(N2+1)*(N3+1)
      CALL SBSRF0(DENS,NTOT,DSURF,LEMPTY)
      IF (LEMPTY) RETURN
C
C Set up matrix for real-space to lattice-index transformation.
C
      XN1=0.99999*FLOAT(N1)
      XN2=0.99999*FLOAT(N2)
      XN3=0.99999*FLOAT(N3)
      DET=BAS(1,1)*BAS(2,2)*BAS(3,3)+BAS(1,2)*BAS(2,3)*BAS(3,1)
     *   +BAS(1,3)*BAS(2,1)*BAS(3,2)-BAS(3,1)*BAS(2,2)*BAS(1,3)
     *   -BAS(3,2)*BAS(2,3)*BAS(1,1)-BAS(3,3)*BAS(2,1)*BAS(1,2)
      IF (ABS(DET).LT.SMALL2) RETURN
      DETNRM=1.0/DET
      MTRX(1,1)=XN1*DETNRM*(BAS(2,2)*BAS(3,3)-BAS(2,3)*BAS(3,2))
      MTRX(1,2)=XN2*DETNRM*(BAS(2,3)*BAS(3,1)-BAS(2,1)*BAS(3,3))
      MTRX(1,3)=XN3*DETNRM*(BAS(2,1)*BAS(3,2)-BAS(2,2)*BAS(3,1))
      MTRX(2,1)=XN1*DETNRM*(BAS(3,2)*BAS(1,3)-BAS(3,3)*BAS(1,2))
      MTRX(2,2)=XN2*DETNRM*(BAS(3,3)*BAS(1,1)-BAS(3,1)*BAS(1,3))
      MTRX(2,3)=XN3*DETNRM*(BAS(3,1)*BAS(1,2)-BAS(3,2)*BAS(1,1))
      MTRX(3,1)=XN1*DETNRM*(BAS(1,2)*BAS(2,3)-BAS(1,3)*BAS(2,2))
      MTRX(3,2)=XN2*DETNRM*(BAS(1,3)*BAS(2,1)-BAS(1,1)*BAS(2,3))
      MTRX(3,3)=XN3*DETNRM*(BAS(1,1)*BAS(2,2)-BAS(1,2)*BAS(2,1))
      CALL SBRCOP(LATICE,ORIG,3)
C
C Some general initialisations.
C
      DDSURF=MAX(ABS(DSURF),SMALL)
      IF (DSURF.LT.0.0) DDSURF=-DDSURF
      GRDSCL(1)=-0.5/(DDSURF*FLOAT(N1))
      GRDSCL(2)=-0.5/(DDSURF*FLOAT(N2))
      GRDSCL(3)=-0.5/(DDSURF*FLOAT(N3))
      COL0=FLOAT(IC1)
      COLSCL=FLOAT(IC2-IC1)
      DO 30 I=1,3
        DXYZ(I,1)=BAS(I,1)/FLOAT(N1)
        DXYZ(I,2)=BAS(I,2)/FLOAT(N2)
        DXYZ(I,3)=BAS(I,3)/FLOAT(N3)
        DDXYZ(I,1,1)=0.0
        DDXYZ(I,1,2)=DXYZ(I,1)
        DDXYZ(I,2,1)=DDXYZ(I,1,1)+DDXYZ(I,1,2)
        DDXYZ(I,2,2)=DXYZ(I,2)
        DDXYZ(I,3,1)=DDXYZ(I,2,1)+DDXYZ(I,2,2)
        DDXYZ(I,3,2)=-DXYZ(I,1)
        DDXYZ(I,4,1)=DDXYZ(I,3,1)+DDXYZ(I,3,2)
        DDXYZ(I,4,2)=-DXYZ(I,2)
        DO 20 J=1,4
          DDXYZ(I,J+4,1)=DDXYZ(I,J,1)
          DDXYZ(I,J+4,2)=DXYZ(I,3)
          DDXYZ(I,J+8,1)=DDXYZ(I,J,1)+DXYZ(I,3)
          DDXYZ(I,J+8,2)=DDXYZ(I,J,2)
  20    CONTINUE
  30  CONTINUE
C
C First paint the edges of the lattice.
C
      DO 60 IFACE=1,3
        I=IFACE
        J=MOD(IFACE,3)+1
        K=MOD(J,3)+1
        IF (IFACE.EQ.1) THEN
          IN=N1
          JN=N2
        ELSEIF (IFACE.EQ.2) THEN
          IN=N2
          JN=N3
        ELSE
          IN=N3
          JN=N1
        ENDIF
        KK=0
        XN=BAS(2,J)*BAS(3,I)-BAS(2,I)*BAS(3,J)
        YN=BAS(3,J)*BAS(1,I)-BAS(3,I)*BAS(1,J)
        ZN=BAS(1,J)*BAS(2,I)-BAS(1,I)*BAS(2,J)
        DNRM=SQRT(XN**2+YN**2+ZN**2+SMALL2)
        PNEYE(1)=DBLE(EYE(1)-0.5*(LATICE(1,I+1)+LATICE(1,J+1)))
        PNEYE(2)=DBLE(EYE(2)-0.5*(LATICE(2,I+1)+LATICE(2,J+1)))
        PNEYE(3)=DBLE(EYE(3)-0.5*(LATICE(3,I+1)+LATICE(3,J+1)))
        DEYE=SNGL(SQRT(PNEYE(1)**2+PNEYE(2)**2+PNEYE(3)**2+DSMAL2))
        XLNORM=DBLE(XN)*PNEYE(1)+DBLE(YN)*PNEYE(2)+DBLE(ZN)*PNEYE(3)
        COSSEE=SNGL(XLNORM)/(DEYE*DNRM)
        IF (COSSEE.LT.0.001) THEN
          KK=N3
          IF (IFACE.EQ.2) KK=N1
          IF (IFACE.EQ.3) KK=N2
          PNEYE(1)=PNEYE(1)+DBLE(BAS(1,K))
          PNEYE(2)=PNEYE(2)+DBLE(BAS(2,K))
          PNEYE(3)=PNEYE(3)+DBLE(BAS(3,K))
          DEYE=SNGL(SQRT(PNEYE(1)**2+PNEYE(2)**2+PNEYE(3)**2+DSMAL2))
          XLNORM=DBLE(XN)*PNEYE(1)+DBLE(YN)*PNEYE(2)+DBLE(ZN)*PNEYE(3)
          COSSEE=-SNGL(XLNORM)/(DEYE*DNRM)
        ENDIF
        IF (COSSEE.GT.0.001) THEN
          XYZ1=FLOAT(KK)*DXYZ(1,K)+LATICE(1,1)
          XYZ2=FLOAT(KK)*DXYZ(2,K)+LATICE(2,1)
          XYZ3=FLOAT(KK)*DXYZ(3,K)+LATICE(3,1)
          DO 50 J1=1,JN
            J0=J1-1
            DO 40 I1=1,IN
              I0=I1-1
              IF (IFACE.EQ.1) THEN
                DLOCAL(1)=DENS(I0,J0,KK)-DSURF
                DLOCAL(2)=DENS(I1,J0,KK)-DSURF
                DLOCAL(3)=DENS(I1,J1,KK)-DSURF
                DLOCAL(4)=DENS(I0,J1,KK)-DSURF
              ELSEIF (IFACE.EQ.2) THEN
                DLOCAL(1)=DENS(KK,I0,J0)-DSURF
                DLOCAL(2)=DENS(KK,I1,J0)-DSURF
                DLOCAL(3)=DENS(KK,I1,J1)-DSURF
                DLOCAL(4)=DENS(KK,I0,J1)-DSURF
              ELSE
                DLOCAL(1)=DENS(J0,KK,I0)-DSURF
                DLOCAL(2)=DENS(J0,KK,I1)-DSURF
                DLOCAL(3)=DENS(J1,KK,I1)-DSURF
                DLOCAL(4)=DENS(J1,KK,I0)-DSURF
              ENDIF
              CALL SBSRF1(DLOCAL,IBSIDE,FRCXYZ)
              IF (IBSIDE.NE.0) THEN
                XYZ(1)=XYZ1+DXYZ(1,I)*FLOAT(I0)
                XYZ(2)=XYZ2+DXYZ(2,I)*FLOAT(I0)
                XYZ(3)=XYZ3+DXYZ(3,I)*FLOAT(I0)
                CALL SBTSF2(XYZ,DXYZ(1,I),DXYZ(1,J),IBSIDE,FRCXYZ,VERT,
     *                      EYE,LIGHT,LSHINE,ITRANS)
              ENDIF
  40        CONTINUE
            XYZ1=XYZ1+DXYZ(1,J)
            XYZ2=XYZ2+DXYZ(2,J)
            XYZ3=XYZ3+DXYZ(3,J)
  50      CONTINUE
        ENDIF
  60  CONTINUE
C
C Step through each "cube" in the lattice, and paint any isosurfaces
C found therein.
C
      X00K=LATICE(1,1)
      Y00K=LATICE(2,1)
      Z00K=LATICE(3,1)
      DO 90 K1=1,N3
        K0=K1-1
        DO 80 J1=1,N2
          J0=J1-1
          X0JK=X00K+DXYZ(1,2)*FLOAT(J0)
          Y0JK=Y00K+DXYZ(2,2)*FLOAT(J0)
          Z0JK=Z00K+DXYZ(3,2)*FLOAT(J0)
          DO 70 I1=1,N1
            I0=I1-1
            DLOCAL(1)=DENS(I0,J0,K0)-DSURF
            DLOCAL(2)=DENS(I1,J0,K0)-DSURF
            DLOCAL(3)=DENS(I1,J1,K0)-DSURF
            DLOCAL(4)=DENS(I0,J1,K0)-DSURF
            DLOCAL(5)=DENS(I0,J0,K1)-DSURF
            DLOCAL(6)=DENS(I1,J0,K1)-DSURF
            DLOCAL(7)=DENS(I1,J1,K1)-DSURF
            DLOCAL(8)=DENS(I0,J1,K1)-DSURF
            CALL SBSRF3(DLOCAL,IVERT,FRCXYZ,ISUMV,ISUMF)
            IF (ISUMV.NE.0) THEN
              XYZ(1)=X0JK+DXYZ(1,1)*FLOAT(I0)
              XYZ(2)=Y0JK+DXYZ(2,1)*FLOAT(I0)
              XYZ(3)=Z0JK+DXYZ(3,1)*FLOAT(I0)
              CALL SBSRF4(DENS,N1,N2,N3,I0,J0,K0,GRDSCL,BAS,GRDCUB)
              CALL SBTSF5(XYZ,DDXYZ,ISUMV,ISUMF,IVERT,FRCXYZ,VERT,EYE,
     *                    LIGHT,LSHINE,ITRANS)
            ENDIF
  70      CONTINUE
  80    CONTINUE
        X00K=X00K+DXYZ(1,3)
        Y00K=Y00K+DXYZ(2,3)
        Z00K=Z00K+DXYZ(3,3)
  90  CONTINUE
      END SUBROUTINE SBTSUR

      SUBROUTINE SBTSF2(XYZ,D1,D2,IB,FRC,VERT,EYE,LIGHT,LSHINE,ITRANS)
C     ----------------------------------------------------------------
C
      REAL    XYZ(*),D1(*),D2(*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      LOGICAL LSHINE
C
      IF (IB.EQ.15) THEN
        DO 10 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=XYZ(I)+D1(I)
          VERT(I,3)=VERT(I,2)+D2(I)
          VERT(I,4)=XYZ(I)+D2(I)
  10    CONTINUE
        CALL SBTSF6(ITRANS,EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.1) THEN
        DO 20 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+FRC(1)*D1(I)
          VERT(I,3)=XYZ(I)+(1.0-FRC(4))*D2(I)
  20    CONTINUE
        CALL SBTSF6(ITRANS,EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.2) THEN
        DO 30 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+FRC(2)*D2(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
  30    CONTINUE
        CALL SBTSF6(ITRANS,EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.4) THEN
        DO 40 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=VERT(I,1)-FRC(3)*D1(I)
          VERT(I,3)=VERT(I,1)-(1.0-FRC(2))*D2(I)
  40    CONTINUE
        CALL SBTSF6(ITRANS,EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.8) THEN
        DO 50 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=VERT(I,1)-FRC(4)*D2(I)
          VERT(I,3)=VERT(I,1)+(1.0-FRC(3))*D1(I)
  50    CONTINUE
        CALL SBTSF6(ITRANS,EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.7) THEN
        DO 60 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+D1(I)
          VERT(I,3)=VERT(I,2)+D2(I)
          VERT(I,4)=VERT(I,3)-FRC(3)*D1(I)
          VERT(I,5)=XYZ(I)+(1.0-FRC(4))*D2(I)
  60    CONTINUE
        CALL SBTSF6(ITRANS,EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.14) THEN
        DO 70 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+D2(I)
          VERT(I,3)=XYZ(I)+D2(I)
          VERT(I,4)=XYZ(I)+(1.0-FRC(4))*D2(I)
          VERT(I,5)=XYZ(I)+FRC(1)*D1(I)
  70    CONTINUE
        CALL SBTSF6(ITRANS,EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.13) THEN
        DO 80 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=XYZ(I)+D2(I)
          VERT(I,3)=XYZ(I)
          VERT(I,4)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,5)=XYZ(I)+D1(I)+FRC(2)*D2(I)
  80    CONTINUE
        CALL SBTSF6(ITRANS,EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.11) THEN
        DO 90 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=XYZ(I)
          VERT(I,3)=XYZ(I)+D1(I)
          VERT(I,4)=VERT(I,3)+FRC(2)*D2(I)
          VERT(I,5)=VERT(I,1)+(1.0-FRC(3))*D1(I)
  90    CONTINUE
        CALL SBTSF6(ITRANS,EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.3) THEN
        DO 100 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=XYZ(I)+D1(I)
          VERT(I,3)=VERT(I,2)+FRC(2)*D2(I)
          VERT(I,4)=XYZ(I)+(1.0-FRC(4))*D2(I)
 100    CONTINUE
        CALL SBTSF6(ITRANS,EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.6) THEN
        DO 110 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+D2(I)
          VERT(I,3)=VERT(I,2)-FRC(3)*D1(I)
          VERT(I,4)=XYZ(I)+FRC(1)*D1(I)
 110    CONTINUE
        CALL SBTSF6(ITRANS,EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.12) THEN
        DO 120 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=XYZ(I)+D2(I)
          VERT(I,3)=VERT(I,2)-FRC(4)*D2(I)
          VERT(I,4)=XYZ(I)+D1(I)+FRC(2)*D2(I)
 120    CONTINUE
        CALL SBTSF6(ITRANS,EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.9) THEN
        DO 130 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=XYZ(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,4)=VERT(I,1)+(1.0-FRC(3))*D1(I)
 130    CONTINUE
        CALL SBTSF6(ITRANS,EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.5) THEN
        DO 140 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+FRC(1)*D1(I)
          VERT(I,3)=XYZ(I)+(1.0-FRC(4))*D2(I)
          VERT(I,4)=XYZ(I)+D1(I)+D2(I)
          VERT(I,5)=VERT(I,4)-FRC(3)*D1(I)
          VERT(I,6)=VERT(I,4)-(1.0-FRC(2))*D2(I)
 140    CONTINUE
        CALL SBTSF6(ITRANS,EYE,3,VERT,LSHINE,LIGHT,0) 
        CALL SBTSF6(ITRANS,EYE,3,VERT(1,4),LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.10) THEN
        DO 150 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+FRC(2)*D2(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,4)=XYZ(I)+D2(I)
          VERT(I,5)=VERT(I,4)-FRC(4)*D2(I)
          VERT(I,6)=VERT(I,4)+(1.0-FRC(3))*D1(I)
 150    CONTINUE
        CALL SBTSF6(ITRANS,EYE,3,VERT,LSHINE,LIGHT,0) 
        CALL SBTSF6(ITRANS,EYE,3,VERT(1,4),LSHINE,LIGHT,0) 
      ENDIF
      END SUBROUTINE SBTSF2

      SUBROUTINE SBTSF5(XYZ,DXYZ,ISV,ISF,IV,FRC,VERT,EYE,LIGHT,LSHINE,
     *                  ITRANS)
C     ----------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV(*),IV4MAP(12)
      LOGICAL LSHINE
      DATA    IV4MAP /12,8,4,3,11,7,6,2,10,9,5,1/
C
      IF (ISV.EQ.1) THEN
        CALL STSF5A(ITRANS,XYZ,DXYZ,FRC,VERT,IV(1),EYE,LSHINE,LIGHT)
      ELSEIF (ISV.EQ.2) THEN
        IF (ISF.EQ.6) THEN
          CALL STSF5A(ITRANS,XYZ,DXYZ,FRC,VERT,IV(1),EYE,LSHINE,LIGHT)
          CALL STSF5A(ITRANS,XYZ,DXYZ,FRC,VERT,IV(2),EYE,LSHINE,LIGHT)
        ELSE
          IJDIF=IV(2)-IV(1)
          IF (IV(1).LE.4) THEN
            IF (IV(2).LE.4) THEN
              K2=IV(1)
              IF (IJDIF.EQ.3) K2=IV(2)
            ELSE
              K2=IV(2)
            ENDIF
          ELSE
            K2=IV(1)+4
            IF (IJDIF.EQ.3) K2=IV(2)+4
          ENDIF
          CALL STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,K2,EYE,LSHINE,LIGHT,1)
        ENDIF
      ELSEIF (ISV.EQ.3) THEN
        IF (ISF.EQ.9) THEN
          DO 10 I=1,3
  10        CALL STSF5A(ITRANS,XYZ,DXYZ,FRC,VERT,IV(I),EYE,LSHINE,LIGHT)
        ELSEIF (ISF.EQ.6) THEN
          DO 20 I1=1,3
            I2=1+MOD(I1,3)
            I=MIN(I1,I2)
            J=MAX(I1,I2)
            K2=0
            IJDIF=IV(J)-IV(I)
            IF (IV(I).LE.4) THEN
              IF (IV(J).LE.4) THEN
                IF (IJDIF.EQ.1) THEN
                  K2=IV(I)
                ELSEIF (IJDIF.EQ.3) THEN
                  K2=IV(J)
                ENDIF
              ELSE
                IF (IJDIF.EQ.4) K2=IV(J)
              ENDIF
            ELSE
              IF (IJDIF.EQ.1) THEN
                K2=IV(I)+4
              ELSEIF (IJDIF.EQ.3) THEN
                K2=IV(J)+4
              ENDIF
            ENDIF
            IF (K2.GT.0) GOTO 1
  20      CONTINUE
   1      CALL STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,K2,EYE,LSHINE,LIGHT,1)
        ELSE
          K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
          CALL STSF5C(ITRANS,XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
        ENDIF
      ELSE
        IF (ISF.EQ.12) THEN
          DO 30 I=1,4
  30        CALL STSF5A(ITRANS,XYZ,DXYZ,FRC,VERT,IV(I),EYE,LSHINE,LIGHT)
        ELSEIF (ISF.EQ.4) THEN
          K4=(IV(1)+IV(2)+IV(3)+IV(4)-6)/4
          IF ((IV(2)-IV(1)).EQ.3) K4=6
          CALL STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT,2)
        ELSEIF (ISF.EQ.6) THEN
          IF (IV(3).LE.4) THEN
            K3=IV(1)+IV(2)+IV(3)-6
            K4=MOD((IV(4)+K3),4)+3*K3
          ELSE
            IF (IV(2).GE.5) THEN
              K3=IV(2)+IV(3)+IV(4)-18
              K4=IV4MAP(MOD((IV(1)+K3),4)+3*K3)
            ELSE
              K4=12+IV(3)-IV(2)
              IF ((IV(1)+IV(2)+IV(3)+IV(4)).EQ.22) K4=29-K4
            ENDIF
          ENDIF
          CALL STSF5D(ITRANS,XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT)
        ELSE
          K4=IV(1)+IV(2)+IV(3)+IV(4)
          IF (K4.EQ.16 .OR. K4.EQ.20) THEN          
            CALL STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,IV(3),EYE,LSHINE,
     *                  LIGHT,1)
            CALL STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,IV(4),EYE,LSHINE,
     *                  LIGHT,1)
          ELSEIF (K4.EQ.18) THEN
            K4A=IV(1)
            IF ((IV(2)-K4A).EQ.3) K4A=4
            K4B=9+MOD(K4A+1,4)
            CALL STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,K4A,EYE,LSHINE,LIGHT,1)
            CALL STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,K4B,EYE,LSHINE,LIGHT,1)
          ELSE
            IF (K4.EQ.14) THEN
              K4A=IV(4)
              K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
            ELSEIF (K4.EQ.22) THEN
              K4A=IV(1)
              K3=IV(2)+IV(3)+IV(4)-5+2*(IV(2)/5+2*(IV(3)/5)+IV(4)/5)
            ELSE
              IF (MOD((IV(1)+IV(2)),2).EQ.0) THEN
                IF (IV(4).EQ.6 .OR. (IV(3)-IV(2)).EQ.2) THEN
                  K4A=IV(2)
                  K3=IV(1)+IV(3)+IV(4)-5+2*(IV(1)/5+2*(IV(3)/5)+IV(4)/5)
                ELSE
                  K4A=IV(1)
                  K3=IV(2)+IV(3)+IV(4)-5+2*(IV(2)/5+2*(IV(3)/5)+IV(4)/5)
                ENDIF
              ELSE
                IF (IV(1).EQ.3 .OR. (IV(3)-IV(2)).EQ.2) THEN
                  K4A=IV(3)
                  K3=IV(1)+IV(2)+IV(4)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(4)/5)
                ELSE
                  K4A=IV(4)
                  K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
                ENDIF
              ENDIF
            ENDIF
            CALL STSF5A(ITRANS,XYZ,DXYZ,FRC,VERT,K4A,EYE,LSHINE,LIGHT)
            CALL STSF5C(ITRANS,XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
          ENDIF
        ENDIF
      ENDIF
      END SUBROUTINE SBTSF5

      SUBROUTINE STSF5A(ITRANS,XYZ,DXYZ,FRC,VERT,IV,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      LOGICAL LSHINE
C
      IF (IV.LE.4) THEN
        J=IV
        K=1+MOD(IV+2,4)
        L=IV+4
      ELSE
        J=IV+4
        K=9+MOD(IV-2,4)
        L=IV
      ENDIF
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
  10  CONTINUE
      CALL SBTSF6(ITRANS,EYE,3,VERT,LSHINE,LIGHT,1) 
      END SUBROUTINE STSF5A

      SUBROUTINE STSF5B(ITRANS,XYZ,DXYZ,FRC,VERT,KK,EYE,LSHINE,LIGHT,LL)
C     ------------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IVL(4,12,2)
      LOGICAL LSHINE
      DATA    IVL /5,6,2,4,6,7,3,1,7,8,4,2,8,5,1,3,9,1,4,12,10,2,1,9,
     *     11,3,2,10,12,4,3,11,12,5,6,10,9,6,7,11,10,7,8,12,11,8,5,9,
     *     5,6,7,8,4,2,10,12,1,3,11,9,4,2,10,12,5,6,7,8,1,3,11,9,24*0/
C
      J=IVL(1,KK,LL)
      K=IVL(2,KK,LL)
      L=IVL(3,KK,LL)
      M=IVL(4,KK,LL)
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
        VERT(I,4)=XYZ(I)+DXYZ(I,M,1)+FRC(M)*DXYZ(I,M,2)
        VERT(I,5)=VERT(I,1)
        VERT(I,6)=0.25*(VERT(I,1)+VERT(I,2)+VERT(I,3)+VERT(I,4))
  10  CONTINUE
      DO 20 I=1,4
        CALL SBRCOP(VERT(1,I),VERT(1,7),6)
        CALL SBTSF6(ITRANS,EYE,3,VERT(1,6),LSHINE,LIGHT,1)
  20  CONTINUE
      END SUBROUTINE STSF5B

      SUBROUTINE STSF5C(ITRANS,XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV3(5,24)
      LOGICAL LSHINE
      DATA    IV3 /5,6,7,3,4,8,5,6,2,3,7,8,5,1,2,6,7,8,4,1,
     *             12,4,2,6,9,4,2,10,9,5,9,1,3,8,12,9,1,3,7,10,
     *             1,3,11,10,6,1,3,11,12,5,4,2,10,11,8,2,4,12,11,7,
     *             10,12,4,1,6,12,10,2,1,5,11,9,1,4,8,1,9,11,7,2,
     *             3,11,9,6,2,3,11,9,5,4,2,10,12,8,3,4,12,10,7,3,
     *             5,6,7,11,12,8,5,6,10,11,7,8,5,9,10,6,7,8,12,9/
C
      J=IV3(1,K3)
      K=IV3(2,K3)
      L=IV3(3,K3)
      M=IV3(4,K3)
      N=IV3(5,K3)
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
        VERT(I,4)=XYZ(I)+DXYZ(I,M,1)+FRC(M)*DXYZ(I,M,2)
        VERT(I,5)=XYZ(I)+DXYZ(I,N,1)+FRC(N)*DXYZ(I,N,2)
        VERT(I,6)=VERT(I,1)
        VERT(I,7)=0.2*(VERT(I,1)+VERT(I,2)+VERT(I,3)+VERT(I,4)+
     *                 VERT(I,5))
  10  CONTINUE
      DO 20 I=1,5
        CALL SBRCOP(VERT(1,I),VERT(1,8),6)
        CALL SBTSF6(ITRANS,EYE,3,VERT(1,7),LSHINE,LIGHT,1)
  20  CONTINUE
      END SUBROUTINE STSF5C

      SUBROUTINE STSF5D(ITRANS,XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV4(6,16)
      LOGICAL LSHINE
      DATA    IV4 /12,9,6,7,3,4,10,9,5,4,3,7,11,10,6,5,4,3,
     *             11,12,5,6,2,3,9,12,8,3,2,6,10,9,5,8,3,2,
     *             10,11,8,5,1,2,12,11,7,2,1,5,9,12,8,7,2,1,
     *             9,10,7,8,4,1,11,10,6,1,4,8,12,11,7,6,1,4,
     *             12,10,6,1,3,8,12,10,7,3,1,5,11,9,6,2,4,8,
     *             11,9,5,4,2,7/
      DATA    VNORM /0.1666666667/
C
      CALL SBRFIL(VERT(1,8),0.0,3)
      DO 20 J=1,6
        K=IV4(J,K4)
        DO 10 I=1,3
          VERT(I,J)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
          VERT(I,8)=VERT(I,8)+VERT(I,J)
  10    CONTINUE
  20  CONTINUE
      CALL SBRCOP(VERT,VERT(1,7),3)
      DO 30 I=1,3
  30    VERT(I,8)=VERT(I,8)*VNORM
      DO 40 I=1,6
        CALL SBRCOP(VERT(1,I),VERT(1,9),6)
        CALL SBTSF6(ITRANS,EYE,3,VERT(1,8),LSHINE,LIGHT,1)
  40  CONTINUE
      END SUBROUTINE STSF5D

      SUBROUTINE SBTSF6(ITRANS,EYE,NV,VERT,LSHINE,LIGHT,INSIDE)
C     ---------------------------------------------------------
C
      REAL             EYE(*),VERT(3,*),LIGHT(*)
      LOGICAL          LSHINE
C
      REAL*8           XLNORM,ZZ,DZZ
      REAL             XW(20),YW(20),MTRX
      LOGICAL          LPS,LCOLOR
      COMMON  /SRFCOM/ GRDCUB(3,8),MTRX(3,3),ORIG(3),XL2,COL0,COLSCL
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C Carry out some initial checks and calculate the coordinates of the 
C projected triangle.
C
      IF (NV.LT.3 .OR. NV.GT.10) RETURN
      SMALL=1.0E-10
      XMIN=+1.0E20
      XMAX=-1.0E20
      YMIN=+1.0E20
      YMAX=-1.0E20
      DO 10 I=1,NV
        CALL SBLIN1(EYE,VERT(1,I),VERT(2,I),VERT(3,I),XW(I),YW(I))
        IF (XW(I).LT.XMIN) THEN
          XMIN=XW(I)
          ILEFT=I
        ENDIF
        IF (YW(I).LT.YMIN) THEN
          YMIN=YW(I)
          JBOTOM=I
        ENDIF
        XMAX=MAX(XW(I),XMAX)
        YMAX=MAX(YW(I),YMAX)
  10  CONTINUE
      IF (XMIN.GE.XTRC .OR. XMAX.LE.XBLC) RETURN
      IF (YMIN.GE.YTRC .OR. YMAX.LE.YBLC) RETURN
C
C Find the outward normal seen by the eye.
C
      AX=VERT(1,2)-VERT(1,1)
      AY=VERT(2,2)-VERT(2,1)
      AZ=VERT(3,2)-VERT(3,1)
      BX=VERT(1,1)-VERT(1,NV)
      BY=VERT(2,1)-VERT(2,NV)
      BZ=VERT(3,1)-VERT(3,NV)
      XN=BY*AZ-AY*BZ
      YN=BZ*AX-AZ*BX
      ZN=BX*AY-AX*BY
      TEN=XN*(EYE(1)-VERT(1,1))+YN*(EYE(2)-VERT(2,1))
     *   +ZN*(EYE(3)-VERT(3,1))
      IF (TEN.LT.0.0) THEN
        XN=-XN
        YN=-YN
        ZN=-ZN
        TEN=-TEN
      ENDIF
C
C Plot the projected triangle.
C
      ITLEVL=MAX(MIN(ITRANS,3),1)
      XLNORM=DBLE(TEN)
      EYENRM=XN*EYE(1)+YN*EYE(2)+ZN*EYE(3)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      DYJ=1.0/DY
      DXI=1.0/DX
      SAFER=0.0001
      IF ((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
        JYMIN=INT((YMIN-YBLC)*DY)+2
        JYMAX=MIN(INT((YMAX-YBLC)*DY)+1,NYP)
        IF (JYMIN.GT.JYMAX) RETURN
        YJ=YBLC+(FLOAT(JYMIN-1)+SAFER)*DYJ
        NVL2=JBOTOM
        NVR2=JBOTOM
        J1=JYMIN
        DO 40 IVERT=1,NV
          IF (YJ.GT.YW(NVL2)) THEN
   1        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVL2)) GOTO 1
            YDIFL=YW(NVL2)-YW(NVL1)
            IF (ABS(YDIFL).LT.SMALL) YDIFL=SMALL
            GRADL=(XW(NVL2)-XW(NVL1))/YDIFL
          ENDIF
          IF (YJ.GT.YW(NVR2)) THEN
   2        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVR2)) GOTO 2
            YDIFR=YW(NVR2)-YW(NVR1)
            IF (ABS(YDIFR).LT.SMALL) YDIFR=SMALL
            GRADR=(XW(NVR2)-XW(NVR1))/YDIFR
          ENDIF
          IF (YW(NVL2).LT.YW(NVR2)) THEN
            J2=MIN(INT((YW(NVL2)-YBLC)*DY)+1,JYMAX)
          ELSE
            J2=MIN(INT((YW(NVR2)-YBLC)*DY)+1,JYMAX)
          ENDIF
          DO 30 J=J1,J2
            IF (J.GE.1) THEN
              JTEST=MOD(J,2)
              IF (ITLEVL.EQ.3 .AND. JTEST.EQ.1) GOTO 29
              XL=XW(NVL1)+GRADL*(YJ-YW(NVL1))
              XR=XW(NVR1)+GRADR*(YJ-YW(NVR1))
              ISTEP=1
              IX1=MAX(INT((XL-XBLC)*DX)+2,1)
              IX2=MIN(INT((XR-XBLC)*DX)+1,NXP)
              IF (IX1.GT.IX2) THEN
                ISTEP=-1
                IX1=MIN(IX1-1,NXP)
                IX2=MAX(IX2+1,1)
              ENDIF
              XI=XBLC+FLOAT(IX1-1)*DXI
              SDXI=FLOAT(ISTEP)*DXI
              DZZ=DBLE(SDXI*XN)
              ZZ=DBLE(EYENRM-XI*XN-YJ*YN)
              K=(J-1)*NXP+IX1
              DO 20 I=IX1,IX2,ISTEP
                ITEST=MOD(I,2)
                IF (ITLEVL.EQ.1) THEN
                  IF ((ITEST+JTEST).EQ.0) GOTO 19
                ELSEIF (ITLEVL.EQ.2) THEN
                  IF ((ITEST+JTEST).EQ.1) GOTO 19
                ELSE
                  IF (ITEST.EQ.1) GOTO 19
                ENDIF
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))
                  IF (INSIDE.EQ.0) THEN
                    GX=XN
                    GY=YN
                    GZ=ZN
                  ELSE
                    CALL SBSF6A(X,Y,Z,ORIG,MTRX,GRDCUB,GX,GY,GZ)
                  ENDIF
                  CALL SBSF6B(EYE,X,Y,Z,GX,GY,GZ,LIGHT,XL2,LSHINE,CLR)
                  SBBUFF(KSTART+K)=COL0+COLSCL*CLR
                ENDIF
  19            XI=XI+SDXI
                ZZ=ZZ-DZZ
                K=K+ISTEP
  20          CONTINUE
            ENDIF
  29        YJ=YJ+DYJ
  30      CONTINUE
          J1=J2+1
          IF (J1.GT.JYMAX) RETURN
  40    CONTINUE
      ELSE
        IXMIN=INT((XMIN-XBLC)*DX)+2
        IXMAX=MIN(INT((XMAX-XBLC)*DX)+1,NXP)
        IF (IXMIN.GT.IXMAX) RETURN
        XI=XBLC+(FLOAT(IXMIN-1)+SAFER)*DXI
        NVL2=ILEFT
        NVR2=ILEFT
        I1=IXMIN
        DO 70 IVERT=1,NV
          IF (XI.GT.XW(NVL2)) THEN
   3        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVL2)) GOTO 3
            XDIFL=XW(NVL2)-XW(NVL1)
            IF (ABS(XDIFL).LT.SMALL) XDIFL=SMALL
            GRADL=(YW(NVL2)-YW(NVL1))/XDIFL
          ENDIF
          IF (XI.GT.XW(NVR2)) THEN
   4        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVR2)) GOTO 4
            XDIFR=XW(NVR2)-XW(NVR1)
            IF (ABS(XDIFR).LT.SMALL) XDIFR=SMALL
            GRADR=(YW(NVR2)-YW(NVR1))/XDIFR
          ENDIF
          IF (XW(NVL2).LT.XW(NVR2)) THEN
            I2=MIN(INT((XW(NVL2)-XBLC)*DX)+1,IXMAX)
          ELSE
            I2=MIN(INT((XW(NVR2)-XBLC)*DX)+1,IXMAX)
          ENDIF
          DO 60 I=I1,I2
            IF (I.GE.1) THEN
              ITEST=MOD(I,2)
              IF (ITLEVL.EQ.3 .AND. ITEST.EQ.1) GOTO 59
              YL=YW(NVL1)+GRADL*(XI-XW(NVL1))
              YR=YW(NVR1)+GRADR*(XI-XW(NVR1))
              ISTEP=1
              JY1=MAX(INT((YL-YBLC)*DY)+2,1)
              JY2=MIN(INT((YR-YBLC)*DY)+1,NYP)
              IF (JY1.GT.JY2) THEN
                ISTEP=-1
                JY1=MIN(JY1-1,NYP)
                JY2=MAX(JY2+1,1)
              ENDIF
              YJ=YBLC+FLOAT(JY1-1)*DYJ
              SDYJ=FLOAT(ISTEP)*DYJ
              DZZ=DBLE(SDYJ*YN)
              ZZ=DBLE(EYENRM-YJ*YN-XI*XN)
              K=(JY1-1)*NXP+I
              KSTEP=ISTEP*NXP
              DO 50 J=JY1,JY2,ISTEP
                JTEST=MOD(J,2)
                IF (ITLEVL.EQ.1) THEN
                  IF ((ITEST+JTEST).EQ.0) GOTO 49
                ELSEIF (ITLEVL.EQ.2) THEN
                  IF ((ITEST+JTEST).EQ.1) GOTO 49
                ELSE
                  IF (JTEST.EQ.1) GOTO 49
                ENDIF
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))
                  IF (INSIDE.EQ.0) THEN
                    GX=XN
                    GY=YN
                    GZ=ZN
                  ELSE
                    CALL SBSF6A(X,Y,Z,ORIG,MTRX,GRDCUB,GX,GY,GZ)
                  ENDIF
                  CALL SBSF6B(EYE,X,Y,Z,GX,GY,GZ,LIGHT,XL2,LSHINE,CLR)
                  SBBUFF(KSTART+K)=COL0+COLSCL*CLR
                ENDIF
  49            YJ=YJ+SDYJ
                ZZ=ZZ-DZZ
                K=K+KSTEP
  50          CONTINUE
            ENDIF
  59        XI=XI+DXI
  60      CONTINUE
          I1=I2+1
          IF (I1.GT.IXMAX) RETURN
  70    CONTINUE
      ENDIF
      END SUBROUTINE SBTSF6

      SUBROUTINE SBFBKG(IC1,IC2,ISHADE)
C     ---------------------------------
C
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Sets the shading for the background. This routine should be
C    called after SBFINT, and COLINT or COLTAB, but before any objects
C    are plotted.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    ISHADE   I*4    I       -      Order of shading (IC1-->IC2 - IC1):
C                                      1 - Bottom to top.
C                                      2 - Left to right.
C                                      3 - Bottom-left to top-right.
C                                      4 - Top-left to bottom-right.
C                                      5 - Bottom, middle and top.
C                                      6 - Left, middle and right.
C                                      7 - Rectangular zoom to centre.
C                                      8 - Elliptical zoom to centre.
C Globals
C    SFTBUF
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBRFIL     Fills a real aray with a constant.
C
C History
C   D. S. Sivia      12 Oct 1995  Initial release.
C-----------------------------------------------------------------------
C
      IF (IBFMOD.EQ.2) RETURN
      NC=IC2-IC1
      NTOT=NXP*NYP
      IF (NC.EQ.0) THEN
        CALL SBRFIL(SBBUFF(KSTART+1),FLOAT(IC1),NTOT)
        RETURN
      ENDIF
      IF (ISHADE.EQ.1) THEN
        COL=FLOAT(IC1)
        DCOL=0.9999*FLOAT(NC)/FLOAT(NYP-1)
        K=KSTART+1
        DO 1 J=1,NYP
           CALL SBRFIL(SBBUFF(K),COL,NXP)
           K=K+NXP
           COL=COL+DCOL
   1    CONTINUE
      ELSEIF (ISHADE.EQ.2) THEN
        COL=FLOAT(IC1)
        DCOL=0.9999*FLOAT(NC)/FLOAT(NXP-1)
        DO 11 I=1,NXP
          DO 10 K=KSTART+I,KSTART+NTOT,NXP
  10        SBBUFF(K)=COL
          COL=COL+DCOL
  11    CONTINUE
      ELSEIF (ISHADE.EQ.3) THEN
        XN=FLOAT(NXP-1)
        YN=FLOAT(NYP-1)
        COL0=FLOAT(IC1)+0.0001*FLOAT(NC)
        DCOL=0.9998*FLOAT(NC)/(XN**2+YN**2)
        K=KSTART+1
        DO 21 J=0,NYP-1
          YNJ=YN*FLOAT(J)
          DO 20 I=0,NXP-1
            SBBUFF(K)=COL0+DCOL*(XN*FLOAT(I)+YNJ)
            K=K+1
  20      CONTINUE
  21    CONTINUE
      ELSEIF (ISHADE.EQ.4) THEN
        XN=FLOAT(NXP-1)
        YN=FLOAT(1-NYP)
        COL0=FLOAT(IC1)+0.0001*FLOAT(NC)
        DCOL=0.9998*FLOAT(NC)/(XN**2+YN**2)
        K=KSTART+1
        DO 31 J=1,NYP
          YNJ=YN*FLOAT(J-NYP)
          DO 30 I=0,NXP-1
            SBBUFF(K)=COL0+DCOL*(XN*FLOAT(I)+YNJ)
            K=K+1
  30      CONTINUE
  31    CONTINUE
      ELSEIF (ISHADE.EQ.5) THEN
        NYP1=1
        NYP2=NYP/2
        COL=FLOAT(IC1)
        DCOL=0.9999*FLOAT(NC)/FLOAT(NYP2-NYP1)
        K=KSTART+1
        DO 41 L=1,2
          IF (L.EQ.2) THEN
            NYP1=NYP2+1
            NYP2=NYP
            COL=FLOAT(IC2)
            DCOL=-0.9999*FLOAT(NC)/FLOAT(NYP2-NYP1)
          ENDIF
          DO 40 J=NYP1,NYP2
             CALL SBRFIL(SBBUFF(K),COL,NXP)
             K=K+NXP
             COL=COL+DCOL
  40       CONTINUE
  41    CONTINUE
      ELSEIF (ISHADE.EQ.6) THEN
        NXP1=1
        NXP2=NXP/2
        COL=FLOAT(IC1)
        DCOL=0.9999*FLOAT(NC)/FLOAT(NXP2-NXP1)
        DO 52 L=1,2
          IF (L.EQ.2) THEN
            NXP1=NXP2+1
            NXP2=NXP
            COL=FLOAT(IC2)
            DCOL=-0.9999*FLOAT(NC)/FLOAT(NXP2-NXP1)
          ENDIF
          DO 51 I=NXP1,NXP2
            DO 50 K=KSTART+I,KSTART+NTOT,NXP
  50          SBBUFF(K)=COL
            COL=COL+DCOL
  51      CONTINUE
  52    CONTINUE
      ELSEIF (ISHADE.EQ.7) THEN
        NXP2=NXP/2+1
        NYP2=NYP/2+1
        XN=1.0/FLOAT(NXP2-1)
        YN=1.0/FLOAT(NYP2-1)
        COL0=FLOAT(IC2)
        DCOL=-0.9999*FLOAT(NC)
        K=KSTART+1
        DO 61 J=1,NYP
          YNJ=ABS(YN*FLOAT(J-NYP2))
          DO 60 I=1,NXP
            XNI=ABS(XN*FLOAT(I-NXP2))
            SBBUFF(K)=COL0+DCOL*MAX(XNI,YNJ)
            K=K+1
  60      CONTINUE
  61    CONTINUE
      ELSEIF (ISHADE.EQ.8) THEN
        NXP2=NXP/2+1
        NYP2=NYP/2+1
        XN=1.0/FLOAT(NXP2-1)
        YN=1.0/FLOAT(NYP2-1)
        COL0=FLOAT(IC2)
        DCOL=-0.9999*FLOAT(NC)
        K=KSTART+1
        DO 71 J=1,NYP
          YNJ=(YN*FLOAT(J-NYP2))**2
          DO 70 I=1,NXP
            XNI=(XN*FLOAT(I-NXP2))**2
            SBBUFF(K)=COL0+DCOL*MIN(XNI+YNJ,1.0)
            K=K+1
  70      CONTINUE
  71    CONTINUE
      ENDIF
      END SUBROUTINE SBFBKG

      SUBROUTINE SBQINF(XLEFT,XRIGHT,YBOT,YTOP,ZBMIN,ZBMAX)
C     -----------------------------------------------------
C
      LOGICAL          LPS,LCOLOR
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Passes back some useful information about the software buffer
C    and canvas to be plotted. All (x,y,z) values are taken to be given
C    in world coordinates.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    XLEFT    R*4    O       -      X-coord of left-hand of window.
C    XRIGHT   R*4    O       -      X-coord of right-hand of window.
C    YBOT     R*4    O       -      Y-coord of bottom of window.
C    YTOP     R*4    O       -      Y-coord of top of window.
C    ZBMIN    R*4    O       -      Minimum Z for distance buffer.
C    ZBMAX    R*4    O       -      Maximum Z for distance buffer.
C
C Globals
C    SFTBUF
C
C External Calls
!




C   SUBROUTINE   DESCRIPTION
C     None.
C
C History
C   D. S. Sivia       6 Jul 1995  Initial release.
C-----------------------------------------------------------------------
C
      XLEFT=XBLC
      XRIGHT=XTRC
      YBOT=YBLC
      YTOP=YTRC
C
      NTOT=NXP*NYP
      ZBMIN=0.0
      ZBMAX=-1.0E20
      DO 10 I=1,NTOT
        ZBMIN=MIN(SBBUFF(I),ZBMIN)
        ZBMAX=MAX(SBBUFF(I),ZBMAX)
  10  CONTINUE
      END SUBROUTINE SBQINF

      SUBROUTINE SBRFIL(X,A,N)
C     ------------------------
C
      REAL X(*)
C
      DO 10 I=1,N
  10    X(I)=A
      END SUBROUTINE SBRFIL

      SUBROUTINE SBRCOP(X,Y,N)
C     ------------------------
C
      REAL X(*),Y(*)
C
      DO 10 I=1,N
  10    Y(I)=X(I)
      END SUBROUTINE SBRCOP
      
      real function zeroin(ax,bx,f,tol)
      real ax,bx,f,tol
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      real  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
c
c  compute eps, the relative machine precision
c
      eps = 1.0
   10 eps = eps/2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0) go to 10
c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (abs(fc) .ge. abs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0*eps*abs(b) + 0.5*tol
      xm = .5*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0) go to 90
c
c is bisection necessary
c
      if (abs(e) .lt. tol1) go to 70
      if (abs(fa) .le. abs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0*xm*s
      q = 1.0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)
c
c adjust signs
c
   60 if (p .gt. 0.0) q = -q
      p = abs(p)
c
c is interpolation acceptable
c
      if ((2.0*p) .ge. (3.0*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/abs(fc))) .gt. 0.0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end function zeroin


