MULTIPLOT START
PAPERFORMAT A4H
*
PLOT   : Absolute Fluxes
*
\LETTERSIZE 0.25
* 
\OFS  2.0 22.0
\INBOX
\NOCOPYRIGHT
\PENDEF=3
\DEFAULTCOLOR=4
\FONT=HELVET
*
\INSTRUCTION EXPORT
\VAR STAR  = BAT99_58

* Name des Models (muss noch ergaenzt werden)
\VAR MODEL = ?

* PHOTOMETRY  (muss noch ergaenzt werden)
\VAR   v   = ?
\VAR   u-b = ?
\VAR   b-v = ?
\VAR   v-r = ?
* IR-PHOTOMETRY from 2MASS (muss noch ergaenzt werden)
\VAR J = ?
\VAR H = ?
\VAR K = ?

* shifts und E_B-V  (muss noch ergaenzt werden)
\VAR shift =  ?
\VAR EBVSMITH = ?

* Pfad fuer das WNL-Gitter mit 40% Wasserstoff
*\VAR PATH = ~praktikum/skripte/d3/models/wnl40/

* Pfad fuer das WNL-Gitter mit 20% Wasserstoff
*\VAR PATH = ~praktikum/skripte/d3/models/wnl20/

* Pfad fuer das WNE-Gitter
*\VAR PATH = ~praktikum/skripte/d3/models/wne/

* UV Beobachtung (muss noch ergaenzt werden)
\VAR IUE-SHORT = ~praktikum/
\VAR IUE-LONG  = ~praktikum/

* optische Beobachtung (ist schon fest verdrahtet)
\VAR OPT       = ~praktikum/skripte/d3/obsdat.dir/Br47.dat
\VAR OPT2      = ~praktikum/skripte/d3/obsdat.dir/Br47.3770-9055

*\VAR vinf  = 1600

* Distance Modulus of the LMC: after HIPPARCOS,
* cf. Madore & Freedman, 1998, ApJ 492, 110
\VAR DM = 18.5
\VAR Mv = ?

\VAR VRAD = 250
\VAR CONVOL = "GAUSS=3.0" 


\EXPR RADFAC = $VRAD / 300000. 
\EXPR RADFAC = 1. - $RADFAC 

\EXPR PATH = $PATH // $MODEL

\EXPR FNAME     = $PATH // /formal.plot
\EXPR FNAMEOUT  = $PATH // /formal.out
\EXPR FNAMEFLUX = $PATH // /wruniq.plot

\EXPR REDD = 1.21 * $EBVSMITH
\EXPR Av = 4.1 * $EBVSMITH

\EXPR NEGREDD = -1 * $REDD
\EXPR DILUTE = 0.4 * $DM

***  Reddening is split into galactic and LMC contribution!
\VAR  REDD_GAL = 0.03
\EXPR REDD = $REDD - $REDD_GAL
\EXPR REDD = $REDD // " LMC"

*
* \VAR-LIST
\INSTRUCTION NO-EXPORT


\IF $v .NE. ?
  \IF $DM .EQ. ?
    \EXPR DM = $v - $Mv
    \EXPR DM = $DM - $Av
    \FORMAT (F5.2) $DM
    \LUN XMAX YMAX -10. U-0.5 0.3  &EM&Tv&M = $Mv mag  
    \NEXTLUN                       &EDM=$DM mag
  \ENDIF

  \IF $Mv .EQ. ?
    \EXPR Mv = $v - $DM
    \EXPR Mv = $Mv - $Av
    \FORMAT (F5.2) $Mv
    \LUN XMAX YMAX -10. U-0.5 0.3  &EDM=$DM mag
    \NEXTLUN                       &EM&Tv&M = $Mv mag  
  \ENDIF                      &EM&Tv&M = $Mv mag  
\ENDIF

* Photometrie-Marken, Magnituden umgerechnet entsprechend der in PRICOLR
* verwendeten Calibration Constants 
* log f-nue = -0.4 * ( 48.64 + mag) 
\BGRLUN COLOR=4
\IF $v .NE. ?
  \CALC logflambda = 0.4 * $v
  \CALC logflambda = -8.404 - $logflambda
  \LUN LOG5160. $logflambda M.0 M.0 0.2 &0$v
\ENDIF

\IF $b-v .NE. ?
  \EXPR b = $b-v + $v
  \FORMAT (F5.2) $b
  \EXPR logflambda = 0.4 * $b
  \EXPR logflambda = -8.240 - $logflambda
  \LUN LOG4270. $logflambda M.0 M.0 0.2 &0$b
\ENDIF

\IF $u-b .NE. ?
  \EXPR u = $u-b + $b
  \FORMAT (F5.2) $u
  \EXPR logflambda = 0.4 * $u
  \EXPR logflambda = -8.103 - $logflambda
  \LUN LOG3650. $logflambda M.0 M.0 0.2 &0$u
\ENDIF

\IF $v-r .NE. ?
  \EXPR r = $v - $v-r
  \FORMAT (F5.2) $r
  \EXPR logflambda = 0.4 * $r
  \EXPR logflambda = -8.535 - $logflambda
  \LUN LOG6000. $logflambda M.0 M.0 0.2 &0$r
\ENDIF

*Photometrie fuer IR nach
*log flamda = -0.4(mag + const.) +log(c/lamda^2) mit c in [A/s]
*CLIGHT=2.9979E18
*und lamda in [A]
*const: in J,H,K,L,M,N   49.488,49.921,50.440,51.605,51.924,53.574
*CENTER FREQUENCIES OF THE COLOR BANDS from PRICOLR
*      1.26E4,1.60E4,2.22E4,3.54E4,4.8E4,1.E5

\IF $J .NE. ?
  \CALC logflambda = -9.519 -  0.4 * $J
  \LUN LOG1.26E4 $logflambda M.0 M.0 0.2 &0$J
\ENDIF

\IF $H .NE. ?
  \CALC logflambda = -9.8998 -  0.4 * $H
  \LUN LOG1.60E4 $logflambda M.0 M.0 0.2 &0$H
\ENDIF

\IF $K .NE. ?
  \CALC logflambda = -10.392 -  0.4 * $K
  \LUN LOG2.22E4 $logflambda M.0 M.0 0.2 &0$K
\ENDIF

\BGRLUN OFF
*
*
*\PEN=4
\LUNINC XMIN YMAX 0. 0.5 0.3 $FNAMEOUT "&E" " MODEL START"
\LUN XMAX YMAX R-0.5 U-0.5 0.6  &E$STAR
\LUN XMAX YMAX -7. U-0.5 0.3  &EE&Tb-v&M=$EBVSMITH
\LUN XMAX YMAX -5. U-1.5 0.3  &Eshift=$shift dex
\IF $VRAD .NE. 0.
\NEXTLUN                    &Ev&Trad&M = $VRAD km/s
\ENDIF
\PEN=1


HEADER:
 X-ACHSE:
 Y-ACHSE:\CENTER\&Elog F&T#l#&M [erg s&H-1&M cm&H-2&M \A&H-1&M]
     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER
 X:  18.0CM      3.05       4.4        .050      .100     0.0000     
 Y:  5.5CM      -17.        -11.5        0.500     1.000     0.0000     
*
* CALCULATION
*
* Continuum flux
*
 N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND X-RANGE 3. 4.4
COMMAND REDDENING $REDD_GAL
COMMAND REDDENING $REDD
COMMAND Y- $DILUTE
COMMAND Y+ $shift
COMMAND XDEX
COMMAND YDEX
COMMAND HLYMANA EBV= $NEGREDD
COMMAND CONVOL GAUSS 8.
COMMAND XLOG
COMMAND YLOG
COMMAND INCLUDE $FNAMEFLUX INCKEY=" PLOT: EMERGENT"
*
* Adding the lines (UV)
 N=?   COLOR=2
COMMAND Y-RANGE 0. 100.
COMMAND XLOG
COMMAND YLOG
COMMAND ARI 2 + 1
COMMAND INCLUDE $FNAME INCKEY="* IUE SHORT"
*
* Adding the lines (OPT)
 N=?   COLOR=2
COMMAND Y-RANGE 0. 100.
COMMAND XLOG
COMMAND YLOG
COMMAND ARI 3 + 1
COMMAND INCLUDE $FNAME INCKEY="* OPT"
*
* OBSERVATION IUE
*
N=? PLOTSYMBOL=5  PEN=1 XYTABLE
COMMAND XLOG
COMMAND Y-RANGE 1.E-16 1.
COMMAND YLOG
COMMAND INCLUDE $IUE-SHORT
*
* OBSERVATION IUE
*
N=? PLOTSYMBOL=5  PEN=1 XYTABLE
COMMAND CONVOL KASTEN=5.
COMMAND XLOG
COMMAND Y-RANGE 1.E-16 1.
COMMAND YLOG
COMMAND INCLUDE $IUE-LONG
*
*
ENDE
*
**********************************************************
*22222222222222222222222222222222222222222222222222222222*
**********************************************************
PLOT   : IUE-SHORT 
*
\LETTERSIZE 0.25
* 
\OFS  2.0 15.5
\INBOX
*\NOCOPYRIGHT
\PENDEF=3
\DEFAULTCOLOR=4
\FONT=HELVET
*
\PEN=1
\LINUN XMIN 1. XMAX 1.  0. 0.
*\LUN  1600. YMIN  0.0 0.4  0.15  &E10x
*\LINUN  1900. 1.  3200. 1.  0. 0.
*\LINUN  1520. .1  1680. .1  0. 0.
*\LINUN  1650. 1.  2000. 1.  0. 0.
*IDENTS
\PEN=2
\INCLUDE ident.dat INCKEY="* IDENT IUE SHORT"
*
*
 HEADER :
 X-ACHSE:
 Y-ACHSE:\CENTER\&EREL. FLUSS 
     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER
 X:  18.0CM     1180.0      1960.0       20.00     100.00     0.0000     
 Y:  5.5CM      0.0000        8.0         0.500     1.000     0.0000     
*
*CALCULATION
*
N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND HLYMANA EBV= $NEGREDD
COMMAND CONVOL GAUSS 8.
COMMAND INCLUDE $FNAME INCKEY="* IUE SHORT"
*
* OBSERVATION
*
N=? PLOTSYMBOL=5  PEN=1 XYTABLE
COMMAND X* $RADFAC
COMMAND INCLUDE $IUE-SHORT
*
* Rektifizierung mit Modell-Kontinuum
* Continuum flux
*
N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND X-RANGE 3. 4.
COMMAND REDDENING $REDD_GAL
COMMAND REDDENING $REDD
COMMAND Y- $DILUTE
COMMAND Y+ $shift
COMMAND XDEX
COMMAND YDEX
COMMAND ARI 2 / 3
COMMAND SKIP
COMMAND INCLUDE $FNAMEFLUX INCKEY=" PLOT: EMERGENT"
*
ENDE
*
***********************************************************************
*333333333333333333333333333333333333333333333333333333333333333333333*
***********************************************************************
*
PLOT   : OPT-BLUE+MID
*
\LETTERSIZE 0.25
* 
\OFS  2.0 9.0
\INBOX
\NOCOPYRIGHT
\PENDEF=3
\DEFAULTCOLOR=4
\FONT=HELVET

\PEN=1
\LINUN  XMIN 1.  XMAX 1.     0. 0.
*
*IDENTS
\PEN=2
\INCLUDE ident.dat INCKEY="* IDENT BLUE"
\INCLUDE ident.dat INCKEY="* IDENT MID"
*
 HEADER :
 X-ACHSE:
 Y-ACHSE:\CENTER\&EREL. FLUSS 
     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER
 X:  18.0CM     3700.0      6000.0       20.00     200.00     0.0000     
 Y:  5.5CM      0.0000        4.0         0.500     1.000     0.0000     

*
*****  OBSERVATION Ollie  *****
*
N=? PLOTSYMBOL=5  PEN=2 XYTABLE
COMMAND X* $RADFAC
COMMAND INCLUDE $OPT
*
****  OBSERVATION CASPEC  ****
*
N=? PLOTSYMBOL=5  PEN=1 XYTABLE
COMMAND SETNAME CASPEC1
COMMAND X-RANGE 5600. 6000.
COMMAND CONVOL $CONVOL
COMMAND X* $RADFAC
COMMAND INCLUDE $OPT2 INCKEY="*ORDER#5"
*
N=? PLOTSYMBOL=5  PEN=1 XYTABLE
COMMAND SETNAME CASPEC2
COMMAND X-RANGE 5600. 6000.
COMMAND CONVOL $CONVOL
COMMAND X* $RADFAC
COMMAND MERGE CASPEC1 CASPEC2
COMMAND SKIP
COMMAND INCLUDE $OPT2 INCKEY="*ORDER#6"

*
*CALCULATION
*
N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND CONVOL $CONVOL
COMMAND INCLUDE $FNAME INCKEY="* OPT"
*

* Inlet 4686
\LUN  4686. YMIN  M0. 0.3  0.15  &E10x

*
*****  OBSERVATION Ollie  *****
*
N=? PLOTSYMBOL=5  PEN=2 XYTABLE
COMMAND X-RANGE 4586. 4786.
COMMAND Y* 0.1
COMMAND X* $RADFAC
COMMAND INCLUDE $OPT
*
*CALCULATION
*
N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND X-RANGE 4586. 4786.
COMMAND Y* 0.1
COMMAND CONVOL $CONVOL
COMMAND INCLUDE $FNAME INCKEY="* OPT"
*

*
ENDE 
*
***********************************************************************
*444444444444444444444444444444444444444444444444444444444444444444444*
***********************************************************************
*
PLOT   : OPT-RED
*
\LETTERSIZE 0.25
* 
\OFS  2.0 2.5
\INBOX
*\NOCOPYRIGHT
\PENDEF=3
\DEFAULTCOLOR=4
\FONT=HELVET
*
\PEN=1
\LINUN  XMIN 1.  XMAX 1.     0. 0.

*IDENTS
\PEN=2
\INCLUDE ident.dat INCKEY="* IDENT RED"
*
*
 HEADER :
 X-ACHSE:
 Y-ACHSE:\CENTER\&EREL. FLUSS 
     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER
 X:   8.0CM     6000.0      7200.0       50.00     250.00     0.0000     
 Y:  5.5CM      0.0000        5.0         0.500     1.000     0.0000     
*
****  OBSERVATION  CASPEC  ****
*
N=? PLOTSYMBOL=5  PEN=2 XYTABLE
COMMAND SETNAME CASPEC1
COMMAND CONVOL $CONVOL
COMMAND X* $RADFAC
COMMAND INCLUDE $OPT2 INCKEY="*ORDER#6"
*
N=? PLOTSYMBOL=5  PEN=1 XYTABLE
COMMAND SETNAME CASPEC2
COMMAND CONVOL $CONVOL
COMMAND X* $RADFAC
COMMAND MERGE CASPEC1 CASPEC2
COMMAND SKIP
COMMAND INCLUDE $OPT2 INCKEY="*ORDER#7"
*
*
*CALCULATION
*
N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND CONVOL $CONVOL
COMMAND INCLUDE $FNAME INCKEY="* OPT"
*
ENDE

**********************************************************
*555555555555555555555555555555555555555555555555555555555
**********************************************************
PLOT   : IUE-LONG 
*
\LETTERSIZE 0.25
* 
\OFS  11.0 2.5
\INBOX
*\NOCOPYRIGHT
\PENDEF=3
\DEFAULTCOLOR=4
\FONT=HELVET
*
\PEN=1
\LINUN XMIN 1. XMAX 1.  0. 0.
*IDENTS
\PEN=2
\INCLUDE ident.dat INCKEY="* IDENT IUE LONG"
*
*
 HEADER :
 X-ACHSE:
 Y-ACHSE:
     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER
 X:  9.0CM     1960.0      3400.        50.00     300.00     0.0000     
 Y:  5.5CM      0.0000        4.0         0.500     1.000     0.0000     
*
*OBSERVATION
*
N=? PLOTSYMBOL=5  PEN=1 XYTABLE
COMMAND CONVOL KASTEN=6.
COMMAND X* $RADFAC
COMMAND INCLUDE $IUE-LONG
*
* Rektifizierung mit Modell-Kontinuum
* Continuum flux
*
N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND X-RANGE 3. 4.
COMMAND REDDENING $REDD_GAL
COMMAND REDDENING $REDD
COMMAND Y- $DILUTE
COMMAND Y+ $shift
COMMAND XDEX
COMMAND YDEX
COMMAND ARI 1 / 2
COMMAND SKIP
COMMAND INCLUDE $FNAMEFLUX INCKEY=" PLOT: EMERGENT"
*
*
*CALCULATION
*
N=?   PLOTSYMBOL=5 SYMBOLSIZE=0.2 COLOR=2
COMMAND CONVOL GAUSS 12.
COMMAND INCLUDE $FNAME INCKEY="* IUE SHORT"
*
ENDE
 

MULTIPLOT END