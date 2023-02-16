PAPERFORMAT A4Q BBROTATE
 PLOT   : Observed Spectrum
\FONT=HELVET
\INSTRUCTION EXPORT
*\INTERACTIVE
\INBOX 
*\PENDEF 3
*\SET_NDATMAX = 400000

**** HERE CHANGE THE NAME OF THE DATASET **************
\VAR file = LWR03570.txt
******************************************************
\VAR ext = $file
\CUTVAR $ext -3:
\ECHO FILE EXTENSION = $ext
\IF $ext .EQ. .txt
\VAR filesh = $file
  \CUTVAR $filesh 1:-4
  \ECHO FILE NAME = $filesh
  \EXPR output = $filesh // .dat
  \EXPR outputnoise = $filesh // .noise
  \EXPR outputbg = $filesh // .bg
\ELSE
  \ECHO FILE NAME = $file
  \EXPR output = $file // .dat
  \EXPR outputnoise = $file // .noise
  \EXPR outputbg = $file // .bg
\ENDIF

\VAR type = $file
\CUTVAR $type 1:2
\IF $type .EQ. SW
  \VAR type = sp
\ELSEIF $type .EQ. LW
  \VAR type = lp
\ENDIF
\ECHO TYPE =  $type
\INSTRUCTION NO-EXPORT

* \INCLUDE $file


 HEADER : $file &4flux &2noise &1background
 X-ACHSE:\CENTER\#l# / \A
 Y-ACHSE:\CENTER\F&T#l#&M [erg/s/cm&H2&M/\A]
     MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER
 X: AUTO   0.       5690.      5703.        10.      100.       0. 
 Y:    0.         -1        20            2           10         0. 

N=? XYTABLE SELECT=-1,1 COLOR=4
COMMAND X- 1
* SHORT-HIGH:
* w(k) = 0.05*k +1150, k = 0, 16499
* Note: put the following COMMAND in the beginning of the COMMAND block:
COMMAND IF $type .EQ. sp
COMMAND X* 0.05
COMMAND X+ 1150.
* LONG-HIGH:
* w(k) = 0.1*k +1900.0, k = 0, 12249
COMMAND ELSEIF $type .EQ. lp
COMMAND X* 0.1
COMMAND X+ 1900.
COMMAND ENDIF
*
COMMAND Y-CUT 0 1
COMMAND WRITE FILE=$output
COMMAND SKIP
COMMAND INCLUDE $file

N=? XYTABLE SELECT=-1,2 COLOR=2
COMMAND X- 1
* SHORT-HIGH:
* w(k) = 0.05*k +1150, k = 0, 16499
* Note: put the following COMMAND in the beginning of the COMMAND block:
COMMAND IF $type .EQ. sp
COMMAND X* 0.05
COMMAND X+ 1150.
* LONG-HIGH:
* w(k) = 0.1*k +1900.0, k = 0, 12249
COMMAND ELSEIF $type .EQ. lp
COMMAND X* 0.1
COMMAND X+ 1900.
COMMAND ENDIF
*
COMMAND Y-CUT 0 0.9
COMMAND WRITE FILE=$outputnoise
COMMAND SKIP
COMMAND INCLUDE $file 

N=? XYTABLE SELECT=-1,3 COLOR=1
COMMAND X- 1
* SHORT-HIGH:
* w(k) = 0.05*k +1150, k = 0, 16499
* Note: put the following COMMAND in the beginning of the COMMAND block:
COMMAND IF $type .EQ. sp
COMMAND X* 0.05
COMMAND X+ 1150.
* LONG-HIGH:
* w(k) = 0.1*k +1900.0, k = 0, 12249
COMMAND ELSEIF $type .EQ. lp
COMMAND X* 0.1
COMMAND X+ 1900.
COMMAND ENDIF
*
COMMAND Y-CUT 0 0.9
COMMAND WRITE FILE=$outputbg
COMMAND SKIP
COMMAND INCLUDE $file 

END

