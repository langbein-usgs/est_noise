#!/bin/sh

#  an extention of EstNoiseAll.sh to do one (or more?) noise models;
#   that noise model is a combination of PLBPx and GM (or GGM) noise

if [ "$#" -lt 1  ]
then
   echo  " "
   echo "  Script determines builds upon EstNoiseAll.sh for a PLBPx + GM noise"
   echo "   Assumes EstNoiseAll has been run successfully"
   echo " Usage:  EstNoiseAdd.sh -d data  "
   exit
fi


#  location of executables
progs=/home/langbein/proglib/est_noise8.00/bin_ifort
#   provide location of GMT -- using either version 5 or 6; 
gmtdir=/usr/local/Cellar/gmt/6.5.0_3/bin

while getopts d: option 
do

     case "$option"
     in
          d)  data=$OPTARG;;
#          t)  type=$OPTARG;;
         \?)  echo "Incomplete set of arguements; type EstNoiseAll.sh without arguments to get documentation"
              exit 1;;
     esac
done

here=`pwd`

cd /tmp/SCRATCH/$data
pwd
#  save this MLE for comparison with BP filtered noise
MLE0x=`grep "PL "  noise_"$data".dat | awk '{print $3}'`
echo $MLE0x

#  test for presence of noise" file

if [ ! -f  noise_"$data".dat ]
then
  echo "the file" noise_"$data".dat " does not exist"
  exit
fi


grep PLBP noise_"$data".dat | sort -g -k 3 | tail -1 > junkPLBP

grep GM noise_"$data".dat | sed '/FOGM/d' | tail -1 > junkGM

cat junkGM
cat junkPLBP

MLE0=`grep RW noise_"$data".dat | head -1 | awk '{print $3}'`
echo MLE0 $MLE0

# parse starting values
wn=`awk '{print $6}' junkGM`
echo $wn
wn100=`awk '{print int(100*$6)}' junkGM`
if [ "$wn100" -lt 1 ]
then
  wn=0.01
fi

plexp1=`awk '{print $7}' junkGM`
plamp1=`awk '{print $8}' junkGM`
gm=`awk '{print $9}' junkGM`
plexp2=`awk '{print $10}' junkGM`
plamp2=`awk '{print $11}' junkGM`
np=`awk '{print $12}' junkPLBP`
bpamp=`awk '{print $13}' junkPLBP`
echo $wn $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp

###########################################################
#####
####    Cycle through four runs of Bandpass filtered noise varying the np (number of poles)
####     Rather than using either FL or RW noise, use PL noise as starting point
####    Select best on based upon MLE and compute its wander
##

###
###  Tabulate the results here for later sorting
rm -f PLBPsort.dat
DiffMLEthres=100     #  sets a threshold for deciding whether the 'fast' version worked

##### 1st model with np=1   #########
bpamp=`echo $bpamp | awk '{print $1+0.25}'`
MOD=PLBP1GM
np=1
file=est2_"$MOD".out
cp est0.in est3.in
cat >> est3.in <<EOF
$wn float
$plamp1 float
$plexp1 float  3.0
$gm float
0.5 2
$np
$bpamp float
$plexp2 fix
$plamp2 fix
0
EOF

"$progs"/est_noise8 < est3.in > $file
cp resid.out resid_"$data"_"$MOD".out
tail -75 $file > est_"$data"_"$MOD".out
cp resid.out resid"$MOD".out
cp max.dat max"$MOD".dat

##  Check to see if BP modeling worked for ModType=f; if not, re-run with ModType=c
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
if [ "$Mtype" = "f" ]
then
  ReDo=0
  MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
  MLEI=`grep "MLE=" $file | tail -1 | awk '{print int($2)}'`
  if [ "$MLE" = "Infinity" -o "$MLE" = "-Infinity" -o "$MLE" = "NaN" -o "$MLE" = "Nan" -o "$MLE" = "nan" -o "$MLE" = "inf -o "$MLEI" -gt 1000000" ]
  then
    ReDo=1
  else 
    echo $MLE $MLE0x
    DiffMLE=`echo $MLE $MLE0x | awk '{print int(sqrt(($1-$2)**2))}'`
    echo $DiffMLE
    if [ "$DiffMLE" -gt "$DiffMLEthres" ]
    then
      ReDo=1
    fi
  fi
  echo ReDo $ReDo
  if [ "$ReDo" -eq 1 ]
  then
   lne=`grep -n "# Noise Model Type" est2.in | sed 's/:/ /' | awk '{print $1 -1 }'`
   sed -n '1,'$lne'p' est2.in > est2a.in
   echo "c  #  Noise Model Type" >> est2a.in
   lne=`expr $lne + 2 `
   sed -n ''$lne',$p' est2.in >> est2a.in
   mv $file "$file"_NoGo
   time "$progs"/est_noise8 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi

MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
AIC=`grep "AIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`
BIC=`grep "BIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`

wh=`grep "white noise=" $file | tail -1 | awk '{printf "%.2f\n", $3}'`
ln=`grep -n "power law noise 1" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
gm=`sed -n ''$ln'p' $file | awk '{printf "%.3f\n", $3}'`
ln=`grep -n "power law noise 2" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp

echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat


dMLE=`echo $MLE $MLE0 | awk '{print $1 - $2}'`

dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> PLBPsort.dat

##### 2nd model with np=2   #########
bpamp=`echo $bpamp | awk '{print $1+0.25}'`
MOD=PLBP2GM
np=2
file=est2_"$MOD".out
cp est0.in est3.in
cat >> est3.in <<EOF
$wn float
$plamp1 float
$plexp1 float  3.0
$gm float
0.5 2
$np
$bpamp float
$plexp2 fix
$plamp2 fix
0
EOF

"$progs"/est_noise8 < est3.in > $file
cp resid.out resid_"$data"_"$MOD".out
tail -75 $file > est_"$data"_"$MOD".out
cp resid.out resid"$MOD".out
cp max.dat max"$MOD".dat

##  Check to see if BP modeling worked for ModType=f; if not, re-run with ModType=c
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
if [ "$Mtype" = "f" ]
then
  ReDo=0
  MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
  MLEI=`grep "MLE=" $file | tail -1 | awk '{print int($2)}'`
  if [ "$MLE" = "Infinity" -o "$MLE" = "-Infinity" -o "$MLE" = "NaN" -o "$MLE" = "Nan" -o "$MLE" = "nan" -o "$MLE" = "inf -o "$MLEI" -gt 1000000" ]
  then
    ReDo=1
  else 
    echo $MLE $MLE0x
    DiffMLE=`echo $MLE $MLE0x | awk '{print int(sqrt(($1-$2)**2))}'`
    echo $DiffMLE
    if [ "$DiffMLE" -gt "$DiffMLEthres" ]
    then
      ReDo=1
    fi
  fi
  echo ReDo $ReDo
  if [ "$ReDo" -eq 1 ]
  then
   lne=`grep -n "# Noise Model Type" est2.in | sed 's/:/ /' | awk '{print $1 -1 }'`
   sed -n '1,'$lne'p' est2.in > est2a.in
   echo "c  #  Noise Model Type" >> est2a.in
   lne=`expr $lne + 2 `
   sed -n ''$lne',$p' est2.in >> est2a.in
   mv $file "$file"_NoGo
   time "$progs"/est_noise8 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi

MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
AIC=`grep "AIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`
BIC=`grep "BIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`

wh=`grep "white noise=" $file | tail -1 | awk '{printf "%.2f\n", $3}'`
ln=`grep -n "power law noise 1" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
gm=`sed -n ''$ln'p' $file | awk '{printf "%.3f\n", $3}'`
ln=`grep -n "power law noise 2" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp

echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat


dMLE=`echo $MLE $MLE0 | awk '{print $1 - $2}'`

dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> PLBPsort.dat

##### 3rd model with np=3   #########
bpamp=`echo $bpamp | awk '{print $1+0.25}'`
MOD=PLBP3GM
np=3
file=est2_"$MOD".out
cp est0.in est3.in
cat >> est3.in <<EOF
$wn float
$plamp1 float
$plexp1 float  3.0
$gm float
0.5 2
$np
$bpamp float
$plexp2 fix
$plamp2 fix
0
EOF

"$progs"/est_noise8 < est3.in > $file
cp resid.out resid_"$data"_"$MOD".out
tail -75 $file > est_"$data"_"$MOD".out
cp resid.out resid"$MOD".out
cp max.dat max"$MOD".dat

##  Check to see if BP modeling worked for ModType=f; if not, re-run with ModType=c
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
if [ "$Mtype" = "f" ]
then
  ReDo=0
  MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
  MLEI=`grep "MLE=" $file | tail -1 | awk '{print int($2)}'`
  if [ "$MLE" = "Infinity" -o "$MLE" = "-Infinity" -o "$MLE" = "NaN" -o "$MLE" = "Nan" -o "$MLE" = "nan" -o "$MLE" = "inf -o "$MLEI" -gt 1000000" ]
  then
    ReDo=1
  else 
    echo $MLE $MLE0x
    DiffMLE=`echo $MLE $MLE0x | awk '{print int(sqrt(($1-$2)**2))}'`
    echo $DiffMLE
    if [ "$DiffMLE" -gt "$DiffMLEthres" ]
    then
      ReDo=1
    fi
  fi
  echo ReDo $ReDo
  if [ "$ReDo" -eq 1 ]
  then
   lne=`grep -n "# Noise Model Type" est2.in | sed 's/:/ /' | awk '{print $1 -1 }'`
   sed -n '1,'$lne'p' est2.in > est2a.in
   echo "c  #  Noise Model Type" >> est2a.in
   lne=`expr $lne + 2 `
   sed -n ''$lne',$p' est2.in >> est2a.in
   mv $file "$file"_NoGo
   time "$progs"/est_noise8 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi

MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
AIC=`grep "AIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`
BIC=`grep "BIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`

wh=`grep "white noise=" $file | tail -1 | awk '{printf "%.2f\n", $3}'`
ln=`grep -n "power law noise 1" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
gm=`sed -n ''$ln'p' $file | awk '{printf "%.3f\n", $3}'`
ln=`grep -n "power law noise 2" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp
 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat


dMLE=`echo $MLE $MLE0 | awk '{print $1 - $2}'`

dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> PLBPsort.dat



##### 4ft model with np=4   #########
bpamp=`echo $bpamp | awk '{print $1+0.25}'`
MOD=PLBP4GM
np=4
file=est2_"$MOD".out
cp est0.in est3.in
cat >> est3.in <<EOF
$wn float
$plamp1 float
$plexp1 float  3.0
$gm float
0.5 2
$np
$bpamp float
$plexp2 fix
$plamp2 fix
0
EOF

"$progs"/oest_noise< est3.in > $file
cp resid.out resid_"$data"_"$MOD".out
tail -75 $file > est_"$data"_"$MOD".out
cp resid.out resid"$MOD".out
cp max.dat max"$MOD".dat

##  Check to see if BP modeling worked for ModType=f; if not, re-run with ModType=c
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
if [ "$Mtype" = "f" ]
then
  ReDo=0
  MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
  MLEI=`grep "MLE=" $file | tail -1 | awk '{print int($2)}'`
  if [ "$MLE" = "Infinity" -o "$MLE" = "-Infinity" -o "$MLE" = "NaN" -o "$MLE" = "Nan" -o "$MLE" = "nan" -o "$MLE" = "inf -o "$MLEI" -gt 1000000" ]
  then
    ReDo=1
  else 
    echo $MLE $MLE0x
    DiffMLE=`echo $MLE $MLE0x | awk '{print int(sqrt(($1-$2)**2))}'`
    echo $DiffMLE
    if [ "$DiffMLE" -gt "$DiffMLEthres" ]
    then
      ReDo=1
    fi
  fi
  echo ReDo $ReDo
  if [ "$ReDo" -eq 1 ]
  then
   lne=`grep -n "# Noise Model Type" est2.in | sed 's/:/ /' | awk '{print $1 -1 }'`
   sed -n '1,'$lne'p' est2.in > est2a.in
   echo "c  #  Noise Model Type" >> est2a.in
   lne=`expr $lne + 2 `
   sed -n ''$lne',$p' est2.in >> est2a.in
   mv $file "$file"_NoGo
   time "$progs"/est_noise8 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi

MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
AIC=`grep "AIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`
BIC=`grep "BIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`

wh=`grep "white noise=" $file | tail -1 | awk '{printf "%.2f\n", $3}'`
ln=`grep -n "power law noise 1" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
gm=`sed -n ''$ln'p' $file | awk '{printf "%.3f\n", $3}'`
ln=`grep -n "power law noise 2" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp

echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat


dMLE=`echo $MLE $MLE0 | awk '{print $1 - $2}'`

dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> PLBPsort.dat

#################
#######  FOR BP noise, get the optimal noise model
###############

MOD=`sort -g -k 2 PLBPsort.dat | tail -1 | awk '{print $1}'`
np=`sort -g -k 2 PLBPsort.dat | tail -1 | awk '{print $9}'`
file=est2_"$MOD".out

#  Get and print the noise model
MLE=`grep "MLE=" $file | tail -1 | awk '{printf "%.2f\n", $2}'`
AIC=`grep "AIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`
BIC=`grep "BIC=" $file | tail -1 | awk '{printf "%.1f\n", $2}'`
wh=`grep "white noise=" $file | tail -1 | awk '{printf "%.2f\n", $3}'`
ln=`grep -n "power law noise 1" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp1=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
gm=`sed -n ''$ln'p' $file | awk '{printf "%.3f\n", $3}'`
ln=`grep -n "power law noise 2" $file | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp2=`sed -n ''$ln'p' $file | awk '{printf "%.2f\n", $2}'`
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.3f\n", $4}'`
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`


cp resid"$MOD".out resid.out
cp max"$MOD".dat max.dat

#  Compute the drift


seed=`cat seed.dat`
if [ "$seed" -eq 0 ]
then
  cat 2945829 > seed.dat
fi
"$progs"/compare_wander8 <<EOF
otr
resid.out
max.dat
201
a
0.5 2
$np
EOF

####   Plot the wander
ls -lt wander.out 
sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/NaN/d' | sed '/Inf/d' > junk
mv junk wander.out
tmin=`head -1 wander.out | awk '{print $1}'`
tmax=`tail -1 wander.out | awk '{print (int($1/1000) + 1)*1000}'`

awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
"$gmtdir"/gmt gmtset MEASURE_UNIT inch
"$gmtdir"/gmt gmtset PAPER_MEDIA letter
"$gmtdir"/gmt gmtset PAPER_MEDIA Custom_612x1000
"$gmtdir"/gmt gmtset HEADER_FONT_SIZE 18p
"$gmtdir"/gmt gmtset HEADER_OFFSET -0.0i
"$gmtdir"/gmt gmtset LABEL_FONT_SIZE 14p
"$gmtdir"/gmt gmtset ANNOT_FONT_SIZE_PRIMARY 12p
"$gmtdir"/gmt gmtset ANNOT_FONT_SIZE_SECONDARY 12p

title=`echo $data  $MOD`
psf=plot2
rm -f plot2.*

awk '{print $1, $2}' wander.out | \
"$gmtdir"/gmt psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:" ":/a2f3:"drift, mm"::."$title":WSen -W3p,255/0/0,solid -Xa1.5 -Ya10.9 -K -P > $psf.ps

awk '{print $1, $3}' wander.out | \
"$gmtdir"/gmt psxy -JX -R -W1p,0/0/0,dashed -K -O -Xa1.5 -Ya10.9 >> $psf.ps
awk '{print $1, $4}' wander.out | \
"$gmtdir"/gmt psxy -JX -R -W1p,0/0/0,dashed -K -O -Xa1.5 -Ya10.9 >> $psf.ps
awk '{print $1, $7}' wander.out | \
"$gmtdir"/gmt psxy -JX -R -W2p,0/0/0,solid -K -O -Xa1.5 -Ya10.9 >> $psf.ps
awk '{print $1, $8}' wander.out | \
"$gmtdir"/gmt psxy -JX -R -W2p,0/0/0,solid -K -O -Xa1.5 -Ya10.9 >> $psf.ps



"$gmtdir"/gmt pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa1.5 -Ya10.9 >> $psf.ps <<EOF
0.08 0.90 8 0 0 ML White Noise $wh
0.06 0.84 8 0 0 ML 1st Power Law
0.08 0.79 8 0 0 ML Index $plexp1
0.08 0.74 8 0 0 ML Amp $plamp1
0.08 0.69 8 0 0 ML GMfreq $gm
0.06 0.63 8 0 0 ML 2nd Power Law
0.08 0.58 8 0 0 ML Index $plexp2
0.08 0.53 8 0 0 ML Amp $plamp2
0.06 0.47 8 0 0 ML BP amp $bpamp
0.08 0.42 8 0 0 ML BP pole $np
EOF

"$gmtdir"/gmt pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa1.5 -Ya10.9 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dMLE
EOF

cp $psf.ps "$here"//"$data"_all2.ps
cp noise_"$data".dat "$here"
cp model_"$data".dat "$here"

echo "composite analysis in " "$here"/"$data"_all2.ps
echo "List of noise models is in " "$here"/noise_"$data".dat
echo "List of time-dependent models is in " "$here"/model_"$data".dat
