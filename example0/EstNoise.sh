#!/bin/sh

#  runs est_noise7.x to determine optimal noise model;
#   then runs compare_wander7 to evaluate the results

if [ "$#" -lt 1  ]
then
   echo  " "
   echo "  Script determines parameters of noise model using maximum likelihood"
   echo " Usage:  EstNoise.sh -d data  -M type "
   echo "  where type is noise model type of the following choices"
   echo "   WH -- white noise only"
   echo "   FL --  Flicker plus white noise"
   echo "   RW --  Random walk plus white noise"
   echo "   PL --  Power law plus white noise"
   echo "   FLRW --  Power law plus white noise"
   echo "   FOGM --  Gauss Markov plus white noise (power law index=2)"
   echo "   GM --   Gauss Markov plus white noise (power law != 2 and estimated"
   echo "   BP1 -- Bandpass filtered noise with 1 pole and RW noise"
   echo "   BP2 -- Bandpass filtered noise with 2 pole and RW noise"
   echo "   BP3 -- Bandpass filtered noise with 3 pole and RW noise"
   echo "   BP4 -- Bandpass filtered noise with 4 pole and RW noise"
   echo "   BP1F -- Bandpass filtered noise with 1 pole and Flicker noise"
   echo "   BP2F -- Bandpass filtered noise with 2 pole and Flicker noise"
   echo "   BP3F -- Bandpass filtered noise with 3 pole and Flicker noise"
   echo "   BP4F -- Bandpass filtered noise with 4 pole and Flicker noise"
   echo " This script requires cleanEst.sh to be run immediately prior to this one"
   echo "   as all of the set-up parameters are done by cleanEst.sh" 
   exit
fi
#  Provide location where the executables are located
progs=/Users/john/proglib/est_noiseBeta/bin

# defaults
nett=otr
netd=otr
tsamDay=1     #  default sampling interval of data --- in days <----  This could be changed if needed
ntype=a
here=`pwd`
echo here $here

while getopts d:M: option 
do

     case "$option"
     in
          d)  data=$OPTARG;;
          M)  modtype=$OPTARG;;
         \?)  echo "Incomplete set of arguements; type EstNoise.sh without arguments to get documentation"
              exit 1;;
     esac
done



cd /tmp/SCRATCH

ntype=`cat modtype.dat`

if [ "$ntype" = "q" ]
then
  ntype=n
fi
nett=`cat nettype.dat`

rm -f est2.in
cp est0.in est2.in
rm -f estn.in
touch estn.in
#  set-up the noise modeling
case "$modtype"
in
     WH ) cat > estn.in <<EOF
1.0  float
0.0001 fix
2    fix
0    fix
.5 2
1
0   fix
2   fix
0   fix
0.0
EOF
       np=1  ;;
     RW ) cat > estn.in <<EOF
1.0  float
1    float
2    fix
0    fix
.5 2
1
0   fix
2   fix
0   fix
0.0
EOF
       np=1  ;;
     FL ) cat > estn.in <<EOF
1.0  float
1    float
1    fix
0    fix
.5 2
1
0   fix
2   fix
0   fix
0.0
EOF
       np=1  ;;
     PL ) cat > estn.in <<EOF
1.0  float
1    float
1    float
0    fix
.5 2
1
0   fix
2   fix
0   fix
0.0
EOF
       np=1  ;;
     FLRW ) cat > estn.in <<EOF
1.0  float
1    float
1    fix
0    fix
.5 2
1
0   fix
2   fix
1   float
0.0
EOF
       np=1  ;;
     FOGM ) cat > estn.in <<EOF
1.0  float
1    float
2    fix
1.0    float
.5 2
1
0   fix
2   fix
0   fix
0.0
EOF
       np=1  ;;
     GM ) cat > estn.in <<EOF
1.0  float
1    float
2    float
10    float
.5 2
1
0   fix
2   fix
0   fix
0.0
EOF
       np=1  ;;
     BP1 ) cat > estn.in <<EOF
1.0  float
1    float
2    fix
0    fix
.5 2
1
0.5   float
2   fix
0   fix
0.0
EOF
       np=1  ;;
     BP2 ) cat > estn.in <<EOF
1.0  float
1    float
2    fix
0    fix
.5 2
2
0.5   float
2   fix
0   fix
0.0
EOF
       np=2  ;;
     BP3 ) cat > estn.in <<EOF
1.0  float
1    float
2    fix
0    fix
.5 2
3
0.5   float
2   fix
0   fix
0.0
EOF
       np=3  ;;
     BP4 ) cat > estn.in <<EOF
1.0  float
1    float
2    fix
0    fix
.5 2
4
0.5   float
2   fix
0   fix
0.0
EOF
       np=4  ;;
     BP1F ) cat > estn.in <<EOF
1.0  float
1    float
1    fix
0    fix
.5 2
1
0.5   float
1   fix
0   fix
0.0
EOF
       np=1  ;;
     BP2F ) cat > estn.in <<EOF
1.0  float
1    float
1    fix
0    fix
.5 2
2
0.5   float
1   fix
0   fix
0.0
EOF
       np=2  ;;
     BP3F ) cat > estn.in <<EOF
1.0  float
1    float
1    fix
0    fix
.5 2
3
0.5   float
2   fix
0   fix
0.0
EOF
       np=3  ;;
     BP4F ) cat > estn.in <<EOF
1.0  float
1    float
1    fix
0    fix
.5 2
4
0.5   float
2   fix
0   fix
0.0
EOF
       np=4  ;;
esac

cat estn.in >> est2.in

cp data.cl data.in

###############
#####
####    Run est_noise7.x to do noise modeling
####
################
echo 772395 > seed.dat
time "$progs"/est_noise7.22 < est2.in > est2.out
cp resid.out resid_"$data"_"$modtype".out
tail -75 est2.out > est_"$data"_"$modtype".out


################
##########
####   Compute the wander/drift and plot
############
############

"$progs"/compare_wander7.01 <<EOF
$netd
resid.out
max.dat
201
$ntype
0.5 2
$np
EOF

####   Plot the wander

tmin=`head -1 wander.out | awk '{print $1}'`
tmax=`tail -1 wander.out | awk '{print (int($1/1000) + 1)*1000}'`

awk '{print $2 }' wander.out > junk
awk '{print $7 }' wander.out >> junk
awk '{print $8 }' wander.out >> junk
min=`sort -g junk | head -1 | awk '{print 0.5*int($1/0.5)}'`
max=`sort -g junk | tail -1 | awk '{print 5*(int($1/5)+1)}'`
if [ "$min" -eq 0 ]
then
  min=`sort -g junk | head -1 |  awk '{print 0.1*int($1/0.1)}'`
fi

echo $tmin/$tmax/$min/$max
gmtset MEASURE_UNIT inch
gmtset PAPER_MEDIA letter
gmtset HEADER_FONT_SIZE 24p
gmtset HEADER_OFFSET -0.1i

title=`echo $data "wander " $modtype`
psf=plot
rm -f plot.ps
echo 
awk '{print $1, $2}' wander.out | \
psxy -JX3.0l -R$tmin/$tmax/$min/$max -Ba1f3:"Period, days":/a2f3:"drift, mm"::."$title":WSen -W3/255/0/0p -X1.5 -Y7 -K -P > $psf.ps

awk '{print $1, $3}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O >> $psf.ps
awk '{print $1, $4}' wander.out | \
psxy -JX -R -W1/0/0/0tap -K -O >> $psf.ps
awk '{print $1, $7}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O >> $psf.ps
awk '{print $1, $8}' wander.out | \
psxy -JX -R -W2/0/0/0p -K -O >> $psf.ps

#  Get and print the noise model
MLE=`grep "MLE=" est2.out | tail -1 | awk '{printf "%.2f\n", $2}'`
AIC=`grep "AIC=" est2.out | tail -1 | awk '{printf "%.1f\n", $2}'`
BIC=`grep "BIC=" est2.out | tail -1 | awk '{printf "%.1f\n", $2}'`
wh=`grep "white noise=" est2.out | tail -1 | awk '{printf "%.2f\n", $3}'`
ln=`grep -n "power law noise 1" est2.out | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp1=`sed -n ''$ln'p' est2.out | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp1=`sed -n ''$ln'p' est2.out | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
gm=`sed -n ''$ln'p' est2.out | awk '{printf "%.3f\n", $3}'`
ln=`grep -n "power law noise 2" est2.out | tail -1 | sed 's/://' | awk '{print $1 + 1}'`
plamp2=`sed -n ''$ln'p' est2.out | awk '{printf "%.2f\n", $2}'`
ln=`expr $ln + 1 `
plexp2=`sed -n ''$ln'p' est2.out | awk '{printf "%.2f\n", $2}'`
bpamp=`grep "Bandpass filter amplitude=" est2.out | tail -1 | awk '{printf "%.2f\n", $4}'`

pstext -JX3.0/3.0 -R0/1/0/1 -K -O >> $psf.ps <<EOF
0.05 0.95 10 0 0 ML MLE $MLE
0.40 0.95 10 0 0 ML AIC $AIC
0.70 0.95 10 0 0 ML BIC $BIC
0.08 0.90 10 0 0 ML White Noise $wh
0.06 0.84 10 0 0 ML 1st Power Law
0.08 0.79 10 0 0 ML Index $plexp1
0.08 0.74 10 0 0 ML Amp $plamp1
0.08 0.69 10 0 0 ML GMfreq $gm
0.06 0.63 10 0 0 ML 2nd Power Law
0.08 0.58 10 0 0 ML Index $plexp2
0.08 0.53 10 0 0 ML Amp $plamp2
0.06 0.47 10 0 0 ML BP amp $bpamp
0.08 0.42 10 0 0 ML BP pole $np
EOF
##########################
###
###   Compute and plot equivalent PSD using the MLE results
###
#################

t1=`head -1 data.in | awk '{printf "%.3f\n", $1 + ($2-1)/365.25}'`
t2=`tail -1 data.in | awk '{printf "%.3f\n", $1 + ($2-1)/365.25}'`
tlen=`echo $t2 $t1 | awk '{printf "%.0f\n", $1 - $2}'`
echo tlen $tlen
"$progs"/psd_calc7 <<EOF
1      #  sampling interval days
16    #  length of time series in year
$plexp1
$plamp1
$gm
$plexp2
$plamp2
$wh
0.5 2
$np
$bpamp
a
EOF


fmin=`head -1 psdcalc.out | awk '{print int(100*$1)/100 -0.01}'`
fmax=200
pmin=`sort -g -k 2 psdcalc.out | head -1 | awk '{print 10*int($2/10)-10}'`
pmax=`sort -g -k 2 psdcalc.out | tail -1 | awk '{print 10*int($2/10)+10}'`
minmax < psdcalc.out
echo $fmin $fmax $pmin $pmax

psxy psdcalc.out -JX2l/2 -R$fmin/$fmax/$pmin/$pmax -Ba1f3:"freq c/yr":/a10f2.5:"mm@+2@+2/c/yr, db":WSen -W3/0/0/0p -K -O -X4.0 >> $psf.ps

#############################################
##########
##########  Plot residuals
#####         only if nett=otr  or nett=gmt
#################

if [ "$nett" = "otr" ]
then
   gmtset PAPER_MEDIA letter
   gmtset HEADER_FONT_SIZE 24p
   gmtset HEADER_OFFSET -0.1i
   title=`echo $data "residual" $modtype`

   R=`awk '{printf "%.3f %.2f\n", $1 +($2-1)/365.25, $3}' resid.out | minmax -I1/5`
   echo $R
   awk '{printf "%.3f %.2f\n", $1 +($2-1)/365.25, $3}' resid.out | \
   psxy -JX6/3 $R -Ba5f1/a10f5:mm::."$title":WSen -Sc0.01 -X-4.4 -Y-5  -O -K >> $psf.ps
fi

if [ "$nett" = "gmt" ]
then
   gmtset PAPER_MEDIA letter
   gmtset HEADER_FONT_SIZE 24p
   gmtset HEADER_OFFSET -0.1i
   title=`echo $data "residual" $modtype`

   tmin=`head -1 resid.out | awk '{print $1}'`
   tmax=`tail -1 resid.out | awk '{print $1}'`
   dmin=`awk '{print $2}' resid.out | sort -g | head -1 | awk '{print 5*int($1/5)-5}'`
   dmax=`awk '{print $2}' resid.out | sort -g | tail -1 | awk '{print 5*int($1/5)+5}'`
   echo $tmin $tmax $dmin $dmax
   awk '{print $1, $2}' resid.out | \
   psxy -JX6T/3 -R"$tmin"/"$tmax"/$dmin/$dmax -Ba5Yf1y/a10f5:mm::."$title":WSen -Sc0.01  -X-4.4 -Y-5 -O -K >> $psf.ps

fi

  mv $psf.ps "$here"/"$data"_wand_"$modtype".ps
  mv est_"$data"_"$modtype".out "$here"/est_"$data"_"$modtype".out
  echo "Results found in"
ls -1 "$here"/"$data"_wand_"$modtype".ps
ls -1 "$here"/est_"$data"_"$modtype".out
          
