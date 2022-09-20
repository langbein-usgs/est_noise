#!/bin/sh

#  runs est_noise7.x to determine optimal noise model;
#   then runs compare_wander7 to evaluate the results

if [ "$#" -lt 1  ]
then
   echo  " "
   echo "  Script determines parameters of noise model using maximum likelihood"
   echo " Usage:  EstNoise.sh -d data  [ -h ] [-S 0/1/2 ] [ -P ] [ -m expmax] [ -B low:high ]"
   echo "   Steps through and evaluates a sequence of noise models likely"
   echo "     to model the background noise of GNSS data"
   echo "   Initially, evaluates random walk and flicker (not the combo)"
   echo "     to establish a null model."
   echo "   Then, proceeds to invoke the power-law, flicker plus random-walk",
   echo "     First order Gauss Markov, and generalized GM noise."
   echo "   For all of these, drift curves are provided and the differences"
   echo "     in Log likelihood relative to the null model are tabulated."
   echo " -S 0/1/2/3  Regulates the types of noise models tested"
   echo "    0 is default and all noise models are tested,"
   echo "      FL, RW, PL, FLRW, two version of Gauss Markov,"
   echo "      and two versions of bandpass filtered noise"
   echo "    1 Tests only FL, RW, PL, FLRW, two version of"
   echo "       Bandpass filtered nose (no Gauss Markov)"
   echo "    2 Tests only FL, RW, PL, FLRW, two version of"
   echo "       Gauss-Markov noise (no Bandpass filtered)"
   echo "    3 Tests only FL, RW, PL, and FLRW noise models"
   echo " -P If P is set, then when Bandpass filtered noise is evaluated"
   echo "     then use PL as instead of FLRW; FLRW is default if -P isn't set"
   echo " -m expmax  Maximum allowed power-law index for first power-law model"
   echo "     default is 3. Used occassionally to limit size of Pl index"
   echo "     Some data appear to have indices > 2.5 which may be unrealistic"
   echo " -B low:high  Sets limits of bandpass filter (: needed to separate the BP filt limits)"
   echo "    default is 0.5:2.0 -- 0.5 c/yr to 2.0 c/yr for modulated, seasonal noise"
   echo " -h --  provides more documentation of algorithm"

   echo " This script requires cleanEst.sh to be run immediately prior to this one"
   echo "   as all of the set-up parameters are done by cleanEst.sh" 
   exit
fi
#  Provide location where the executables are located
#source /home/langbein/.bashrc
progs=/Users/john/proglib/est_noise20151128/bin
progs=/Users/john/Desktop/est_noise20151217/bin
progs=/Users/john/proglib/est_noise20160201/bin
progs=/Users/john/proglib/est_noiseBeta/bin
#progs=/home/langbein/proglib/est_noiseBetaX/bin
# defaults
nett=otr
netd=otr
help=no
bpother=no   #  By default, use FLRW to as the base noise model to test the presence of Bandpass filtered noise
ntype=a
tsamDay=0     #  default sampling interval of data --- in days <----  This could be changed if needed
Nnoise=0      # Default on noise models being tested; Does 8 models
expmax=3    #  default value of power law index

bpfilt=0
skip=0  #  This is for testing; normaly set to zero

nsim=201    #  Number of simulations for compare_wander.


###


while getopts d:hS:Pm:B: option 
do

     case "$option"
     in
          d)  data=$OPTARG;;
          h)  help=yes;;
          S) Nnoise=$OPTARG;;
          P) bpother=yes;;
          m) expmax=$OPTARG;;
          B) bpfilt=$OPTARG;;
         \?)  echo "Incomplete set of arguements; type EstNoiseAll.sh without arguments to get documentation"
              exit 1;;
     esac
done

if [ "$help" = "yes" ]
then
  rm -f /tmp/junk
  cat > /tmp/junk <<EOF

	 More discussion of EstNoiseAll.sh

This script provides a model which can be modified to measure both
the background noise and parameters that describe a time dependent
function of a GNSS time series. It is assumed that the data are
measure once per day.  In particular, EstNoiseAll.sh demonstrates
the interrelationship between est_noise7.xx and compare_wander7.

Prior to running EstNoiseAll+ one needs to run cleanEst.sh which
does preliminary fitting of prescribed functions to input time
series. It assumes that the background noise is equal parts white,
flicker, and random walk. One should play with the various available
functions (offsets, trends), and editing.  Once a "reasonable" set
of functions is described, cleanEst.sh will remove outliers using
the bust_5 program. Importantly, cleanEst.sh provides a clean data
set in /tmpSCRATCH/$data/data.cl and the initial half of the set-up to
drive est_noise7 (/tmp/SCRATCH/$data/est0.in).  est0.in contains the
configuration used by est_noise7 to model the time dependence.
EstNoiseAll+ will append to est0.in the necessary configuration for
the noise modeling.

EstNoiseAll.sh evaluates 8 different noise models allowed in
est_noise7. These are flicker noise (FL), random-walk noise (RW),
power-law noise (PL), a combination of FL and RW nose,  First Order
Gauss Markov noise (FOGM), Generalized Gauss Markov noise (GM),
bandpass filter noise plus either FL or RW noise (BPx), and bandpass
filtered noise plus FLRW noise (FLRWBPx) or plus PL (PLBPx), where
"x" is the number of poles used in the bandpass filter. The bandpass
filter spans the seasonal (1-year) noise sources.

Initially, the script finds the optimum FL noise with est_noise7
and its corresponding drift spectrum using compare_wander7. The
values of the white noise and flicker amplitude are saved with a
'fl' suffix.

Second: The optimal RW noise is determined using est_noise7 and its
corresponding drift using compare_wander7. The values of the white
noise and RW amplitude are saved with a 'rw' suffix.

Third: The values of log-likelihood, MLE, for both FL and RW are
compared and the largest one is determined to be the "better" noise
model and identified as the "null" noise model for comparisons with
the more complex noise models.  The values MLE, AIC, and BIC are
computed relative to those of the so-called null model and those
values are placed on the two drift plots for comparison; Those with
zero are the "null" results.

Forth: Using the values of noise for the null model, these are used
a the initial guess to a power-law (PL) noise model.  est_noise7
is run to estimate the terms of the PL noise; compare_wander7 then
is used to compute the drift spectrum and that is plotted. The
parameters of the PL noise are saved with a 'pl' suffix for later
use with bandpass filtered noise.  Also, the  MLE, AIC and BIC
relative to the null model are computed and placed on the drift
plot.

Fifth: Likewise, evaluate the combination of flicker and random-walk
noise (FLRW).

Sixth:  Evaluate the First Order Gauss-Markov noise (FOGM).  Use
the values of noise for the RW model as an initial guess.  Note
that FOGM uses a power-law index of 2 (as with RW). The additional
parameter is the GM frequency in radians/yr.

Seventh: Evaluate the generalize Gauss-Markov noise (GM).  Use the
values of noise for the PL model as an initial guess.

Eighth: Evaluate bandpass filtered noise (centered at 1 c/yr) for
either a power law index of 1 (FL) or 2 (RW) as specified by null
model.  This is done by using an initial guess of the null-model,
for 4 different runs of est_noise7. Each run of est_noise7 uses a
different number of bandpass filter poles (np), from 1 to 4.  At
the conclusion of the 4 runs, the result with the highest MLE is
selected, and that model is evaluated with compare_wander7.  Note
that both AIC and BIC are re-evaluated as est_noise7 does not account
for the number of poles being a free parameter.

Lastly,  Evaluate bandpass filtered noise (centered at 1 c/yr)
superimposed on either FLRW or PL noise.  This is done by using an
initial guess of the PL-model or FLRW-model, for 4 different runs
of est_noise7. Each run of est_noise7 uses a different number of
bandpass filter poles (np), from 1 to 4.  At the conclusion of the
4 runs, the result with the highest MLE is selected, and that model
is evaluated with compare_wander7.  Note that both AIC and BIC are
re-evaluated as est_noise7 does not account for the number of poles
being a free parameter.  Note that the choice to use either FLRW,
the default, or PL is toggled with the -P option.


With BP filtered noise and the use of the "fast" option of additive
noise, it is possible that est_noise7 will not successfully converge.
The script will sense that problem and re-start the BP noise modeling but
using the slower, cholesky decomposition method to invert the data covariance.

By default, the script runs through all eight types of noise models.
Often, neither Gauss Markov or Bandpass filter noise is needed and
the script provides an option not to evaluate these models of noise.

At the conclusion, there will be 8 plots of date-drift superimposed
on the statistical ranges of simulated data having the same noise
as the estimated from the actual data.  The red line represents the
drift of the data adjusted for its time-dependence.  The dashed
black curve represents the range of 68% of the simulations. Likewise,
the solid black curve is the range for 95% of the simulations. In
principle, the noise model with the most positive MLE could be the
best.  However, MLE does not take into account the number of unknowns.
These are taken into accounted with AIC and BIC in slightly different
forms. The noise model with the most negative going AIC and BIC
could be designated as "best".  However, one should look closely
at the drift of the data relative to the simulations; if the data
drift exceeds the ranges (especially grossly) of the simulations,
then it could suggest that something is not modeled correctly. One
should examine the fit of the time-dependent model to the data
(obtained from cleanEst.sh) for potential problems.

	

EOF
  cat /tmp/junk
  exit

fi


here=`pwd`
nett=$netd


###  All work happens in /tmp/SCRATCH/$data

cd /tmp/SCRATCH/"$data"
echo 772395 > seed.dat
ntype=`cat modtype.dat`
#bpfilt=`cat bpfilt.dat`
if [ "$bpfilt"  = "0" ] 
then
  echo 0.5 2.0 > bpfilt.dat
else
  echo $bpfilt | sed 's/:/ /' > bpfilt.dat
fi
bpfilt=`cat bpfilt.dat`


if [ "$ntype" = "q" ]
then
  ntype=n
fi
nett=`cat nettype.dat`
netd=`cat nettype.dat`

if [ "$Nnoise" -gt 3 ]
then
  echo "Invalid argument for -S, must be 0, 1, 2, or 3"
  echo "Bailing"
  exit
fi
#   This is the cleaned data created by cleanEst.sh

cp data.cl data.in

####  Create the inputs needed for the various noise models
#   There are 8 of them!!

cat > est_RW.in <<EOF
-99999.  float
-99999.    float
2    fix  $expmax
0    fix
$bpfilt
1
0   fix
2   fix
0   fix
0.0
EOF
cat > est_FL.in <<EOF
-99999.  float
-99999.    float
1    fix  $expmax
0    fix
$bpfilt
1
0   fix
2   fix
0   fix
0.0
EOF
cat > est_PL.in <<EOF
sig  float
plamp1    float
plexp1    float  $expmax
0    fix
$bpfilt
1
0   fix
2   fix
0   fix
0.0
EOF
cat > est_FLRW.in <<EOF
sig  float
plamp1    float
1    fix  $expmax
0    fix
$bpfilt
1
0   fix
2   fix
plamp2   float
0.0
EOF

cat > est_FOGM.in <<EOF
sig  float
plamp1    float
2    fix  $expmax
0.2    float
$bpfilt
1
0   fix
2   fix
0   fix
0.0
EOF

cat > est_GM.in <<EOF
sig  float
plamp1    float
plexp1    float  $expmax
0.2    float
$bpfilt
1
0   fix
2   fix
0   fix
0.0
EOF

cat > est_BP.in <<EOF
sig  float
plamp1    float
plexp1    fix  $expmax
0   fix
$bpfilt
np
0.5   float
2   fix
0   fix
0.0
EOF

cat > est_BPPL.in <<EOF
sig  float
plamp1    float
plexp1    float  $expmax
0 fix
$bpfilt
np
0.5   float
2   fix
0   fix
0.0
EOF

cat > est_BPFLRW.in <<EOF
sig  float
plamp1    float
1    fix  $expmax
0 fix
$bpfilt
np
0.5   float
2   fix
plamp2   float
0.0
EOF

#
## Listing of all noise models
rm -f noise_"$data".dat
cat > noise_"$data".dat << EOF
Model  MLE  AIC BIC  sig plexp1 plamp1 GMrad plexp2 plamp2 poles BPamp
EOF
rm -f model_"$data".dat
cat > model_"$data".dat << EOF
Model	rate	intercept	365	182	offset
EOF

DiffMLEthres=100     #  sets a threshold for deciding whether the 'fast' version worked

####   Establish null noise model by running both Flicker and Random walk
###
###  Flicker

MOD=FL
np=1
rm -f est2.in
cp est0.in est2.in
cat est_"$MOD".in >> est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]      
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi

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
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat


#  Compute the drift

cp resid"$MOD".out resid.out
cp max"$MOD".dat max.dat
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

####   Plot the wander
ls -lt wander.out 
sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' | sed '/NaN/d'  > junk
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
gmtset MEASURE_UNIT inch
gmtset PAPER_MEDIA letter
gmtset PAPER_MEDIA Custom_612x1000
gmtset HEADER_FONT_SIZE 18p
gmtset HEADER_OFFSET -0.2i
gmtset LABEL_FONT_SIZE 14p
gmtset ANNOT_FONT_SIZE_PRIMARY 12p
gmtset ANNOT_FONT_SIZE_SECONDARY 12p

title=`echo $data  $MOD`
psf=plot
rm -f plot.ps
echo 
awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:" ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa1.5 -Ya10.9 -K -P > $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya10.9 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya10.9 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya10.9 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya10.9 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa1.5 -Ya10.9 >> $psf.ps <<EOF
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
MLEFL=$MLE

AICFL=$AIC
BICFL=$BIC
whfl=$wh
plampfl=$plamp1
plexpfl=$plexp1
##0.05 0.95 8 0 0 ML MLE $MLE
##0.40 0.95 8 0 0 ML AIC $AIC
##0.70 0.95 8 0 0 ML BIC $BIC


###########   Do Random Walk   ##########
MOD=RW
np=1
rm -f est2.in
cp est0.in est2.in
cat est_"$MOD".in >> est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi

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
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
MLERW=$MLE
AICRW=$AIC
BICRW=$BIC
whrw=$wh
plamprw=$plamp1
plexprw=$plexp1
########################################
#
#
#  Establish whether FL or RW is null
#
####################

MLERW100=`echo $MLERW | awk '{print int(100*$1)}'`
MLEFL100=`echo $MLEFL | awk '{print int(100*$1)}'`
echo RW stats: $MLERW $AICRW $BICRW $whrw $plamprw $plexprw
echo FL stats: $MLEFL $AICFL $BICFL $whrw $plampfl $plexpfl
echo RW $MLERW100 $MLERW $AICRW $BICRW $whrw $plamprw $plexprw > junk
echo FL $MLEFL100 $MLEFL $AICFL $BICFL $whrw $plampfl $plexpfl >> junk
sort -g -k 2 junk
null=`sort -g -k 2 junk | tail -1 | awk '{print $1}'`
echo NULL is $null
if [ "$null" = "FL" ]
then
  MLE0=$MLEFL
  AIC0=$AICFL
  BIC0=$BICFL
  wh0=$whfl
  plamp0=$plampfl
  plexp0=$plexpfl
else
  MLE0=$MLERW
  AIC0=$AICRW
  BIC0=$BICRW
  wh0=$whrw
  plamp0=$plamprw
  plexp0=$plexprw
fi
echo NULL $null $MLE0 $AIC0 $BIC0 $wh0 $plamp0 $plexp0

dmle=`echo $MLEFL $MLE0 | awk '{print $1-$2}'`
daic=`echo $AICFL $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BICFL $BIC0 | awk '{print $1-$2}'`
pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa1.5 -Ya10.9 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
EOF



dmle=`echo $MLERW $MLE0 | awk '{print $1-$2}'`
daic=`echo $AICRW $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BICRW $BIC0 | awk '{print $1-$2}'`

cp resid"$MOD".out resid.out
cp max"$MOD".dat max.dat
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF
sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:" ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa1.5 -Ya8.0 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya8.0 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya8.0 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya8.0 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya8.0 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa1.5 -Ya8.0 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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

#######  Done with null noise model   ##########

#########  NEXT IS PL model    ###############
MOD=PL
np=1
rm -f est2.in
cp est0.in est2.in
cat est_"$MOD".in >> est2.in
#  Use noise parameters from null model as the initial guess
sed 's/sig/'$wh0'/' est2.in | sed 's/plamp1/'$plamp0'/' | sed 's/plexp1/'$plexp0'/' > junk.in
mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`
whpl=$wh
plamppl=$plamp1
plexppl=$plexp1
#  save this MLE for comparison with BP filtered noise
MLE0x=$MLE

cp resid"$MOD".out resid.out
cp max"$MOD".dat max.dat
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:" ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa1.5 -Ya5.1 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya5.1 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya5.1 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya5.1 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya5.1 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa1.5 -Ya5.1 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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


####################
###   Do FLRW Noise   

##############

MOD=FLRW
np=1
rm -f est2.in
cp est0.in est2.in
cat est_"$MOD".in >> est2.in
#  Use the values of the null noise model as the initial guess.
if [ "$null" = "FL" ]
then
  sed 's/sig/'$whfl'/' est2.in | sed 's/plamp1/'$plampfl'/' | sed 's/plexp1/'$plexpfl'/' | sed 's/plexp2/2.0/' | sed 's/plamp2/1.0/' > junk.in
else
  sed 's/sig/'$whrw'/' est2.in | sed 's/plamp2/'$plamprw'/' | sed 's/plexp2/'$plexprw'/' | sed 's/plexp1/1.0/' | sed 's/plamp1/1.0/' > junk.in
fi
mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi

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
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`
whflrw=$wh
plamp1flrw=$plamp1
plamp2flrw=$plamp2

cp resid"$MOD".out resid.out
cp max"$MOD".dat max.dat
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:"Period, days ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa1.5 -Ya2.2 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya2.2 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa1.5 -Ya2.2 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya2.2 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa1.5 -Ya2.2 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa1.5 -Ya2.2 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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

if [ "$Nnoise" -eq 0 -o "$Nnoise" -eq 2 ]
then

#########  NEXT IS FOGM model    ###############
MOD=FOGM
np=1
rm -f est2.in
cp est0.in est2.in
cat est_"$MOD".in >> est2.in
# use randow-walk noise model as the initial guess
sed 's/sig/'$whrw'/' est2.in | sed 's/plamp1/'$plamprw'/'  > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi

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
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

cp resid"$MOD".out resid.out
cp max"$MOD".dat max.dat
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:" ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa5.0 -Ya10.9 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya10.9 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya10.9 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya10.9 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya10.9 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa5.0 -Ya10.9 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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



#########  NEXT IS Generalized GM model    ###############
MOD=GM
np=1
rm -f est2.in
cp est0.in est2.in
cat est_"$MOD".in >> est2.in
#  Use the results of the PL noise as an initial guess
sed 's/sig/'$whpl'/' est2.in | sed 's/plamp1/'$plamppl'/' | sed 's/plexp1/'$plexppl'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
bpamp=`grep "Bandpass filter amplitude=" $file | tail -1 | awk '{printf "%.2f\n", $4}'`
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

cp resid"$MOD".out resid.out
cp max"$MOD".dat max.dat
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:" ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa5.0 -Ya8.0 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya8.0 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya8.0 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya8.0 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya8.0 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa5.0 -Ya8.0 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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

fi   #   End testing/doing Nnoise for GM

if [ "$Nnoise" -eq 0 -o "$Nnoise" -eq 1 ]
then

###########################################################
#####
####    Cycle through four runs of Bandpass filtered noise varying the np (number of poles)
####    First, though select whether FL or RW noise is used as the basis (Note, later, we'll do PL noise
####    Select best on based upon MLE and compute its wander
##

###
###  Tabulate the results here for later sorting
rm -f BPsort.dat

MOD=BP1
np=1
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BP.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$wh0'/' est2.in | sed 's/plamp1/'$plamp0'/' | sed 's/plexp1/'$plexp0'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
    echo $MLE $MLE0
    DiffMLE=`echo $MLE $MLE0 | awk '{print int(sqrt(($1-$2)**2))}'`
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi

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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> BPsort.dat


#######################################
MOD=BP2
np=2
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BP.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$wh0'/' est2.in | sed 's/plamp1/'$plamp0'/' | sed 's/plexp1/'$plexp0'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
    echo $MLE $MLE0
    DiffMLE=`echo $MLE $MLE0 | awk '{print int(sqrt(($1-$2)**2))}'`
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi


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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`


echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> BPsort.dat


#######################################
MOD=BP3
np=3
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BP.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$wh0'/' est2.in | sed 's/plamp1/'$plamp0'/' | sed 's/plexp1/'$plexp0'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
    echo $MLE $MLE0
    DiffMLE=`echo $MLE $MLE0 | awk '{print int(sqrt(($1-$2)**2))}'`
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat



  fi
fi

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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`


echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> BPsort.dat

#######################################
MOD=BP4
np=4
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BP.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$wh0'/' est2.in | sed 's/plamp1/'$plamp0'/' | sed 's/plexp1/'$plexp0'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
    echo $MLE $MLE0
    DiffMLE=`echo $MLE $MLE0 | awk '{print int(sqrt(($1-$2)**2))}'`
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`


echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> BPsort.dat

######################################
#######  FOR BP noise, get the optimal noise model
##################

MOD=`sort -g -k 2 BPsort.dat | tail -1 | awk '{print $1}'`
np=`sort -g -k 2 BPsort.dat | tail -1 | awk '{print $9}'`
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
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:" ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa5.0 -Ya5.1 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya5.1 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya5.1 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya5.1 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya5.1 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa5.0 -Ya5.1 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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


cat BPsort.dat

if [ "$bpother" = "yes" ]
then
###########################################################
#####
####    Cycle through four runs of Bandpass filtered noise varying the np (number of poles)
####     Rather than using either FL or RW noise, use PL noise as starting point
####    Select best on based upon MLE and compute its wander
##

###
###  Tabulate the results here for later sorting
rm -f PLBPsort.dat



MOD=PLBP1
np=1
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPPL.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whpl'/' est2.in | sed 's/plamp1/'$plamppl'/' | sed 's/plexp1/'$plexppl'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> PLBPsort.dat

MOD=PLBP2
np=2
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPPL.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whpl'/' est2.in | sed 's/plamp1/'$plamppl'/' | sed 's/plexp1/'$plexppl'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> PLBPsort.dat

MOD=PLBP3
np=3
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPPL.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whpl'/' est2.in | sed 's/plamp1/'$plamppl'/' | sed 's/plexp1/'$plexppl'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> PLBPsort.dat

MOD=PLBP4
np=4
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPPL.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whpl'/' est2.in | sed 's/plamp1/'$plamppl'/' | sed 's/plexp1/'$plexppl'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
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
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:"period, days ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa5.0 -Ya2.2 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya2.2 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya2.2 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya2.2 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya2.2 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa5.0 -Ya2.2 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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


else    ###  Rather than BPxPL, do BPxFLRW

###########################################################
#####
####    Cycle through four runs of Bandpass filtered noise varying the np (number of poles)
####     Rather than using either FL or RW noise, use FLRW noise as starting point
####    Select best on based upon MLE and compute its wander
##

###
###  Tabulate the results here for later sorting
rm -f FLRWBPsort.dat

MOD=FLRWBP1
np=1
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPFLRW.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whflrw'/' est2.in | sed 's/plamp1/'$plamp1flrw'/' | sed 's/plamp2/'$plamp2flrw'/' > junk.in

mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> FLRWBPsort.dat

MOD=FLRWBP2
np=2
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPFLRW.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whflrw'/' est2.in | sed 's/plamp1/'$plamp1flrw'/' | sed 's/plamp2/'$plamp2flrw'/' > junk.in


mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> FLRWBPsort.dat

MOD=FLRWBP3
np=3
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPFLRW.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whflrw'/' est2.in | sed 's/plamp1/'$plamp1flrw'/' | sed 's/plamp2/'$plamp2flrw'/' > junk.in


mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> FLRWBPsort.dat

MOD=FLRWBP4
np=4
rm -f est2.in
cp est0.in est2.in
sed 's/np/'$np'/' est_BPFLRW.in > est_"$MOD".in
cat est_"$MOD".in >> est2.in

sed 's/sig/'$whflrw'/' est2.in | sed 's/plamp1/'$plamp1flrw'/' | sed 's/plamp2/'$plamp2flrw'/' > junk.in


mv junk.in est2.in
echo $MOD
file=est2_"$MOD".out
if [ "$skip" -eq 0 ]
then
  time "$progs"/est_noise7.30 < est2.in > $file
  cp resid.out resid_"$data"_"$MOD".out
  tail -75 $file > est_"$data"_"$MOD".out
  cp resid.out resid"$MOD".out
  cp max.dat max"$MOD".dat
fi
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
   time "$progs"/est_noise7.30 < est2a.in > $file
   cp resid.out resid_"$data"_"$MOD".out
   tail -75 $file > est_"$data"_"$MOD".out
   cp resid.out resid"$MOD".out
   cp max.dat max"$MOD".dat
  fi
fi
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
Mtype=`grep "ModType" $file | grep "calling funmin" | awk '{print $2}'`
echo $MOD $Mtype $MLE $AIC $BIC $wh $plexp1 $plamp1 $gm $plexp2 $plamp2 $np $bpamp >> noise_"$data".dat 
echo $MOD $Mtype > zz
paste zz model.dat >> model_"$data".dat
dmle=`echo $MLE $MLE0 | awk '{print $1-$2}'`
daic=`echo $AIC $AIC0 | awk '{print $1-$2}'`
dbic=`echo $BIC $BIC0 | awk '{print $1-$2}'`

##  Modify both AIC and BIC to account for extra noise parameter, np, not accounted for
###   in est_noise
AIC=`echo $AIC | awk '{printf "%.2f\n", $1 + 2 }'`
npts=`wc -l resid_"$data"_"$MOD".out | awk '{print $1}'`
BIC=`echo $BIC $npts | awk '{printf "%.2f\n", $1 + log($2)}'`

echo $MOD $MLE $AIC $BIC $wh $plamp1 $plexp1 $bpamp $np >> FLRWBPsort.dat


#################
#######  FOR BP noise, get the optimal noise model
###############

MOD=`sort -g -k 2 FLRWBPsort.dat | tail -1 | awk '{print $1}'`
np=`sort -g -k 2 FLRWBPsort.dat | tail -1 | awk '{print $9}'`
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
"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
$nsim
$ntype
$bpfilt
$np
EOF

sed '/Nan/d' wander.out | sed '/nan/d' | sed '/NAN/d' | sed '/Inf/d' > junk
mv junk wander.out
awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
#min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
#min10=`sort -g junk | head -1 | awk '{print 5*int($1/0.5)}'`
min10=`sort -g junk | head -1 | awk '{printf "%.0f\n", (10*log($1)/log(10))- 1 }'`
min=`echo $min10 | awk '{print 10**($1/10)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
#if [ "$min10" -eq 0 ]
#then
#  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
#fi

echo $tmin/$tmax/$min/$max
title=`echo $data  $MOD`

awk '{print $1, $2}' wander.out | \
psxy -JX2.5l/2.0l -R$tmin/$tmax/$min/$max -Ba1f3:"period, days ":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -Xa5.0 -Ya2.2 -K -O >> $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya2.2 >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O -Xa5.0 -Ya2.2 >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya2.2 >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O -Xa5.0 -Ya2.2 >> $psf.ps



pstext -JX2.5/2.0 -R0/1/0/1 -K -O -Xa5.0 -Ya2.2 >> $psf.ps <<EOF
0.05 0.95 8 0 0 ML @~d@~MLE $dmle
0.40 0.95 8 0 0 ML @~d@~AIC $daic
0.70 0.95 8 0 0 ML @~d@~BIC $dbic
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

fi   #   Finsh choice between doing BPxPL or PBxFLRW

fi   #  end testing Nnoise  for BP filtered


cp plot.ps "$here"/"$data"_all.ps
cp noise_"$data".dat "$here"
cp model_"$data".dat "$here"

echo "composite analysis in " "$here"/"$data"_all.ps
echo "List of noise models is in " "$here"/noise_"$data".dat
echo "List of time-dependent models is in " "$here"/model_"$data".dat
exit


