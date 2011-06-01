#!MC 1200
# Created by Tecplot 360 build 12.2.0.9077
$!VarSet |MFBD| = '/home/jgraham/Data/lesgo/cubes-ibm/case6_2'
$!READDATASET  '"|MFBD|/mu__0102.dat" "|MFBD|/mu__0202.dat" "|MFBD|/mu__0302.dat" "|MFBD|/mu__0402.dat" "|MFBD|/mu__0502.dat" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = XYLINE
  VARNAMELIST = '"V1" "V2"'
$!ALTERDATA 
  EQUATION = '{H}=15'
$!ALTERDATA 
  EQUATION = '{Ur}=5.04985'
$!ALTERDATA 
  EQUATION = '{y}=V1/{H}'
$!ALTERDATA 
  EQUATION = '{<u>}=V2'
$!ALTERDATA 
  EQUATION = '{Us} = 0.4'
$!ALTERDATA 
  EQUATION = '{alpha_U} = {Us}/{Ur}'
$!ALTERDATA 
  EQUATION = '{<u>s}=0.0'
$!ALTERDATA  [1]
  EQUATION = '{<u>s} = {alpha_U} * {<u>} - 0.3  '
$!ALTERDATA  [2]
  EQUATION = '{<u>s} = {alpha_U} * {<u>} + 0.3  '
$!ALTERDATA  [3]
  EQUATION = '{<u>s} = {alpha_U} * {<u>} + 1.3  '
$!ALTERDATA  [4]
  EQUATION = '{<u>s} = {alpha_U} * {<u>} + 1.7  '
$!ALTERDATA  [5]
  EQUATION = '{<u>s} = {alpha_U} * {<u>} + 2.3  '
$!WRITEDATASET  "|MFBD|/u_scaled_z-0.5.dat"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  INCLUDEAUTOGENFACENEIGHBORS = YES
  ASSOCIATELAYOUTWITHDATAFILE = NO
  VARPOSITIONLIST =  [3-9]
  BINARY = NO
  USEPOINTFORMAT = NO
  PRECISION = 9
  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT
$!RemoveVar |MFBD|
