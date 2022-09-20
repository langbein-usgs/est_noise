#!/bin/sh

#   A script that removes outliers from GPS (or similar) data
#   At first cut, just requires a time series of data
#     Addition information might include 
#      times of offsets
#      intervals of rate changes
#      periods of sinusoid beyand the 365 and 182 day periods
#      exponential and/or Omori law function
#    And known periods of data to be deleted
#
#    Program uses est_noise7 to estimate trends and other functions
#     from raw data
#    Then uses program adjust_1 to remove those trends and other functions
#    Feeds the detrended data into bust_5 which removes outliers after
#     high-pass filtering with a running median
#   Then uses adjust_1 to add back the trends and other functions
#
#   At various states, data are plotted using the gmt package.
#     AT THIS TIME (Dec 2015), plotting is done if data formated using 'otr'
#
#   Written by John Langbein, USGS, Menlo Park, CA
#

if [ "$#" -lt 1  ]
then
   echo  " "
   echo "  Script removes outliers from GPS time series"
   echo " Usage:  cleanEst.sh -d data -f format_type [-b start_time -b end_time -O off_file -E edit_file -T exp_file -S sine_period_file -n noiseType -i sample_interval_days -M noise_model_file -n a_or_c_f_n]"
   echo "   options -b and -e not implimented "
   echo "   only accepts otr and gmt formats"
   echo "   noiseType is a additive (default) or q, quadrature"
   echo "   exp file can accept a list of rate change intervals, and/or"
   echo "     exponentials/omori law trends"
   echo "   default sample interval is 1 day"
   echo "   -M (file) option to use noise model that differs from the default"
   echo "     of a unit of each of white, flicker and random walk"
   echo "     that is a nominal characterization of GPS noise"
   echo "     for daily observations"
   echo "   -S (file) list additional periods (other than 365.25 and 182.625 days)"
   echo "      to estimate their amplitudes/phases"
   echo "   -n default is a -- additive covariance that will use f or c"
   echo "      f -  fast covariance; c - cholesky decomposition or"
   echo "      n -- legacy mode using quadrature covariance"
   echo "      the default is a "
   exit
fi

#  Provide location where the executables are located
#source /home/langbein/.bashrc
progs=/Users/john/proglib/est_noise20151128/bin
progs=/Users/john/Desktop/est_noise20151217/bin
progs=/Users/john/proglib/est_noise20160201/bin
progs=/Users/john/proglib/est_noiseBeta/bin
#progs=/home/langbein/proglib/est_noiseBeta/bin
#  defaults 


tstart=
tstop=
offfile=NO
trdfile=NO
edtfile=NO
nett=otr
netd=otr
ntype=a
tsamDay=1     #  default sampling interval of data --- in days <----  This could be changed if needed
bpfilt=`echo 0.5 2`   #  set-up limits of bp filtered noise
noisefile=NO
sineFile=NO

while getopts d:f:E:O:T:n:i:M:S: option 
do

     case "$option"
     in
          d)  data=$OPTARG;;
          f)  netd=$OPTARG;;
          E)  edtfile=$OPTARG;;
          O)  offfile=$OPTARG;;
          T)  trdfile=$OPTARG;;
          n)  ntype=$OPTARG;;
          i)  tsamDay=$OPTARG;;
          M)  noisefile=$OPTARG;;
          S)  sineFile=$OPTARG;;
         \?)  echo "Incomplete set of arguements; type cleanEst.sh without arguments to get documentation"
              exit 1;;
     esac
done
here=`pwd`
mkdir /tmp/SCRATCH
mkdir /tmp/SCRATCH/"$data"

nett=$netd
if [ "$ntype" = "q" ]
then
  ntype=n
fi

############
##
##   Initial data preparation : copy data and remove known problematic data if there is an edit file
#
###############
cp $data /tmp/SCRATCH/"$data"/data.in


cd /tmp/SCRATCH/"$data"
#   Store the noise model type for future use
rm -f modtype.dat
echo $ntype > modtype.dat
rm -f nettype.dat
echo $netd > nettype.dat

echo $bpfilt > bpfilt.dat

if [ "$edtfile" != "NO" ]
then
#  echo will edit data

  if [ ! -s "$here"/"$edtfile" ]
  then
     echo $edtfile "  Does not exist"
     ls -l "$here"/"$edtfile"
     echo "Bailing!"
     exit
  fi
  npts=`wc -l data.in | awk '{print $1}'`

  if [ "$nett" = "otr" ]
  then
    nedt=`wc -l "$here"/"$edtfile" | awk '{print $1}'`
    n=1
    while [ "$n" -le "$nedt" ]
    do
      t1=`sed -n ''$n'p' "$here"/"$edtfile" | awk '{print $1*1000+$2}'`
      t2=`sed -n ''$n'p' "$here"/"$edtfile" | awk '{print $3*1000+$4}'`
      rm -f tmp
      awk '$1*1000+$2 < '$t1' {print $1, $2, $3, $4}' data.in > tmp
      awk '$1*1000+$2 > '$t2' {print $1, $2, $3, $4}' data.in >> tmp
      mv tmp data.in
      n=`expr $n + 1`
    done
  fi


  if [ "$nett" = "gmt" ]
  then
    nedt=`wc -l "$here"/"$edtfile" | awk '{print $1}'`
    n=1
    awk '{print $1}' data.in | sed 's/://g' | sed 's/-//g' | sed 's/T//g' > time
    paste time data.in > tmp1

    while [ "$n" -le "$nedt" ]
    do
      t1=`sed -n ''$n'p' "$here"/"$edtfile" | awk '{print $1}' | sed 's/://g' | sed 's/-//g' | sed 's/T//g'`
      t2=`sed -n ''$n'p' "$here"/"$edtfile" | awk '{print $2}'| sed 's/://g' | sed 's/-//g' | sed 's/T//g'`
      echo $t1 $t2

 
      awk '$1 < '$t1' {print $1, $2, $3, $4}' tmp1 > tmp
      awk '$1 > '$t2' {print $1, $2, $3, $4}' tmp1 >> tmp
      mv tmp tmp1
      n=`expr $n + 1 `
    done
    awk '{print $2, $3, $4}' tmp1 > data.in
  fi

  npts1=`wc -l data.in | awk '{print $1}'`
  ndel=`echo $npts $npts1 | awk '{print $1 -$2}'`
  echo "Number of point deleted " $ndel
fi



######################################################
###
####   This section organizes the rather long input needed to drive est_noise7.x   #####
#
#    It will get information needed to characterize the time-dependent functions used to
#       model the data
#
##########################

rm -f est0.in

if [ -z "$tstart" ]
then
   if [ "$nett" = "otr" ]
   then
     tstart=`head -1 data.in | awk '{print $1, $2 }' `
   fi
   if [ "$nett" = "gmt" ]
   then
     tstart=`head -1 data.in | awk '{print $1 }' `
   fi
fi
echo $tstart
if [ -z "$tstop" ]
then
   if [ "$nett" = "otr" ]
   then
     tstop=`tail -1 data.in | awk '{print $1, $2 }' `
   fi
   if [ "$nett" = "gmt" ]
   then
     tstop=`tail -1 data.in | awk '{print $1}' `
   fi
fi
echo $tstop

cat > est0.in <<EOF
$nett
1     # number of baselines
$tstart $tstop        # interval to analyze
y
EOF

#   get rate change information
rm -f rate.in tmp
if [ "$trdfile" != "NO" ]
then
  echo trend file exists
  if [ ! -s "$here"/"$trdfile" ]
  then
     echo $trdfile "  Does not exist"
     ls -l "$here"/"$trdfile"
     echo "Bailing!"
     exit
  fi
  nrate=0
  ntot=`wc -l "$here"/"$trdfile" | awk '{print $1}'`
  n=1
  while [ "$n" -le "$ntot" ]
  do
    tp=`sed -n ''$n'p' "$here"/"$trdfile" | awk '{print $1}'`
    if [ "$tp" = "r"  ]
    then
      if [ "$nett" = "otr" ]
      then
        sed -n ''$n'p' "$here"/"$trdfile" | awk '{print $2, $3, $4, $5}' >> tmp
      fi
      nrate=`expr $nrate + 1 `
    fi
    n=`expr $n + 1 `
  done
  echo $nrate "# number of rate changes" > rate.in
  cat tmp >> rate.in
else  
  echo "0   # number of rate changes" > rate.in
fi


cat rate.in >> est0.in

#  Get sinusoids
rm -f sea.in
touch sea.in
cat > sea.in <<EOF
365.25
182.625
EOF

if [ "$sineFile" != "NO"  ]
then
   awk '{print $1}'  "$here"/"$sineFile" >> sea.in
fi

wc -l sea.in | awk '{print $1 }' > tmp
mv sea.in tmp1
cat tmp tmp1 > sea.in
rm tmp*
cat sea.in >> est0.in

#   Get offsets

rm -f off.in
if [ "$offfile" != "NO" ]
then
  if [ ! -s "$here"/"$offfile" ]
  then
     echo $offfile "  Does not exist"
     ls -l "$here"/"$offfile"
     echo "Bailing!"
     exit
  fi

  wc -l "$here"/"$offfile" | awk '{print $1}' > off.in
  cat "$here"/"$offfile" >> off.in
else

  echo "0   #  Number of offsets" > off.in

fi
cat off.in >> est0.in


#  Get exponential/Omori
rm -f exp.in tmp
if [ "$trdfile" != "NO" ]
then
  echo trend file exists
  if [ ! -s "$here"/"$trdfile" ]
  then
     echo $trdfile "  Does not exist"
     ls -l "$here"/"$trdfile"
     echo "Bailing!"
     exit
  fi
  nexp=0
  ntot=`wc -l "$here"/"$trdfile" | awk '{print $1}'`
  n=1
  while [ "$n" -le "$ntot" ]
  do
    tp=`sed -n ''$n'p' "$here"/"$trdfile" | awk '{print $1}'`
    if [ "$tp" = "e" -o "$tp" = "m" ]
    then
      if [ "$nett" = "otr" ]
      then
        texp=`sed -n ''$n'p' "$here"/"$trdfile" | awk '{print $2, $3}'`
        initial=`sed -n ''$n'p' "$here"/"$trdfile" | awk '{print $4, $5}'`
      fi
      cat >> tmp <<EOF
$texp
$initial
$tp
EOF
      nexp=`expr $nexp + 1`
    fi

    n=`expr $n + 1`
  done
  echo $nexp > exp.in
  cat tmp >> exp.in
else
  echo "0    # Number of exponentials " > exp.in
fi

cat exp.in >> est0.in

#  List the data
cat >> est0.in <<EOF
$netd
data.in
0     #  Number of auxilary files
$tsamDay      #  sampling interval
$ntype    # Noise Model Type
n
0     #  decimate type
EOF

###   This completes the time-dependent modeling section to drive est_noise
#
###  Next, add the noise modeling ---  For data clean-up, I will assume
#      that noise is equal components white, flicker and random walk
#
######
cp est0.in est1.in

if [ "$noisefile" != "NO" ]
then
#  query noise file
   wn=`grep "wn" "$here"/$noisefile | awk '{print $2}'`
   plamp1=`grep "plamp1" "$here"/$noisefile | awk '{print $2}'`
   plexp1=`grep "plexp1" "$here"/$noisefile | awk '{print $2}'`
   gm=`grep "gm" "$here"/$noisefile | awk '{print $2}'`
   np=`grep "np" "$here"/$noisefile | awk '{print $2}'`
   bpamp=`grep "bpamp" "$here"/$noisefile | awk '{print $2}'`
   plamp2=`grep "plamp2" "$here"/$noisefile | awk '{print $2}'`
   plexp2=`grep "plexp2" "$here"/$noisefile | awk '{print $2}'`
   cat >> est1.in <<EOF
$wn  fix    #white noise
$plamp1   fix    # plamp1
$plexp1   fix  3.0  # plexp1
$gm   fix    # GM 
0.50000        2.00000    #  bandpass filter elements
$np    # number of poles
$bpamp fix    #  BP amplitude
$plexp2   fix   #plexp2
$plamp2   fix   #plam2
0.0   # additive white noise
EOF
else

#  assumes combo of white, flicker and RW
  cat >> est1.in <<EOF
1.0   fix    #white noise
1.0   fix    # plamp1
1.0   fix   3.0  # plexp1
0.0   fix    # GM 
0.50000        2.00000    #  bandpass filter elements
1    # number of poles
 0.000 fix    #  BP amplitude
2.0   fix   #plexp2
1.0   fix   #plam2
0.0   # additive white noise
EOF

fi

###############
#####
####    Run est_noise7.x to detrend the data
####
################
echo 7723957 > seed.dat
time "$progs"/est_noise7.30 < est1.in > est1.out
tr -cd '\11\12\15\40-\176' < est1.out > est1a.out
mv est1a.out est1.out

cp resid.out resid1_"$data".out
#tail -75 est1.out > est_"$data".out


##############################
######
######   Remove trends from data
#####      First, tabulate results from est_noise and put these into adj.in
####       Second, run the adjust_1 program
########
#######################

rm -f adj.in

Rate=`grep "Rate in units per year" est1.out | tail -1 | awk '{print $6}'`
cat > adj.in <<EOF
R
$Rate
EOF


nsea=`head -1 sea.in`
n=1
while [ "$n" -le "$nsea" ]
do
   ln=`expr $n + 1`
   per=`sed -n ''$ln'p' sea.in`
   echo $per
   grep "Period of" est1.out | tail -"$nsea" | sed -n ''$n'p' > junk
   camp=`awk '{print $7}' junk`
   samp=`awk '{print $12}' junk`
   echo $camp $samp
   cat >> adj.in <<EOF
s
$tstart
$per $camp $samp
EOF
   n=$ln
done


noff=`head -1 off.in`
n=1
while [ "$n" -le "$noff" ]
do
   ln=`expr $n + 1`
   toff=`sed -n ''$ln'p' off.in`
   if [ "$netd" = "otr" ]
   then
     off=`grep "Offset number" est1.out | tail -"$noff" | sed -n ''$n'p' | awk '{print $8}'`
   fi
   if [ "$netd" = "gmt" ]
   then
     off=`grep "Offset number" est1.out | tail -"$noff" | sed -n ''$n'p' | awk '{print $7}'`
   fi
   echo OFFSET $toff $off
   cat >> adj.in <<EOF
o
$toff
$off
EOF
   n=$ln
done


nexp=`head -1 exp.in| awk '{print $1}'`

n=1
while [ "$n" -le "$nexp" ]
do
  grep "Exponential number" est1.out | tail -"$nexp" |  sed -n ''$n'p'
  ln=`echo $n | awk '{print 1 + $1*3}'`
  tp=`sed -n ''$ln'p' exp.in`
  ln=`echo $n | awk '{print 1 + $1*3-2}'`
  texp=`sed -n ''$ln'p' exp.in`
  if [ "$nett" = "otr" ]
  then
    amp=`grep "Exponential number" est1.out | tail -"$nexp" |  sed -n ''$n'p' | awk '{print $8}'`
    tau=`grep "Exponential number" est1.out | tail -"$nexp" |  sed -n ''$n'p' | awk '{print $14}'`
  fi
  cat >> adj.in <<EOF
$tp
$texp
$tau $amp
EOF
  n=`expr $n + 1 `
done

nrate=`head -1 rate.in | awk '{print $1}'`

n=1
while [ "$n" -le "$nrate" ]
do
   rate=`grep "Rate change number" est1.out | tail -"$nrate" | sed -n ''$n'p' | awk '{print $12}'`
   ln=`expr $n + 1 `
   trate=`sed -n ''$ln'p' rate.in`
   cat >> adj.in <<EOF
r
$trate
$rate
EOF
   n=`expr $n + 1`
done
cat adj.in


######
###
##   Run the adjust_1 program
#
##########

"$progs"/adjust_1 <<EOF
$netd
data.in
$nett
adj.in
d
EOF

#############################################
##########
##########  Plot the detrended data
#####         only if nett=otr  or nett=gmt
#################

if [ "$nett" = "otr" ]
then
   gmtset MEASURE_UNIT inch
   gmtset PAPER_MEDIA letter
   gmtset HEADER_FONT_SIZE 24p
   gmtset HEADER_OFFSET -0.1i
   title=`echo $data "detrended"`
   psf=plot
   rm -f plot.ps
   R=`awk '{printf "%.3f %.2f\n", $1 +($2-1)/365.25, $3}' data.in | minmax -I1/5`
   echo $R
   awk '{printf "%.3f %.2f\n", $1 +($2-1)/365.25, $3}' data.in | \
   psxy -JX6/3 $R -Ba5f1/a10f5:mm::."$title":WSen -Sc0.03 -X1.1 -Y7 -P -K > $psf.ps
fi

if [ "$nett" = "gmt" ]
then
   gmtset MEASURE_UNIT inch
   gmtset PAPER_MEDIA letter
   gmtset HEADER_FONT_SIZE 24p
   gmtset HEADER_OFFSET -0.1i
   title=`echo $data "detrended"`
   psf=plot
   rm -f plot.ps
   tmin=`head -1 data.in | awk '{print $1}'`
   tmax=`tail -1 data.in | awk '{print $1}'`
   dmin=`awk '{print $2}' data.in | sort -g | head -1 | awk '{print 5*int($1/5)-5}'`
   dmax=`awk '{print $2}' data.in | sort -g | tail -1 | awk '{print 5*int($1/5)+5}'`
   echo $tmin $tmax $dmin $dmax
   awk '{print $1, $2}' data.in | \
   psxy -JX6T/3 -R"$tmin"/"$tmax"/$dmin/$dmax -Ba5Yf1y/a10f5:mm::."$title":WSen -Sc0.03 -X1.1 -Y7 -P -K > $psf.ps
fi



###############################
###########
####
#####   Remove outliers 
###
###################  

"$progs"/bust_5 <<EOF
$netd
data.in
data.cl
$tsamDay
90     #   running median filter time
4      #   IQR for rejections
EOF

cp data.cl data.detrend.cl
###################
#######
#######
#####       Add back the trends for further analysis
#######
#############
"$progs"/adjust_1 <<EOF
$netd
data.cl
$nett
adj.in
a
EOF

#####
###
#   re-run est_noise to get a better model (after outliers have been remove) for plots of residuals

cp data.cl data.in

time "$progs"/est_noise7.30 < est1.in > est1.out
tr -cd '\11\12\15\40-\176' < est1.out > est1a.out
mv est1a.out est1.out

cp resid.out resid2_"$data".out
tail -75 est1.out > est_"$data".out


#   replot data with outliers removed
if [ "$nett" = "otr" ]
then

   R=`awk '{printf "%.3f %.2f\n", $1 +($2-1)/365.25, $3}' resid.out | minmax -I1/5`
#   echo $R
   title=`echo $data "detrended, cleaned"`
   awk '{printf "%.3f %.2f\n", $1 +($2-1)/365.25, $3}' resid.out | \
   psxy -JX $R -Ba5f1/a10f5:mm::."$title":WSen -Sc0.03 -X0 -Y-4.0  -O -K >> $psf.ps
fi

#   replot data with outliers removed
if [ "$nett" = "gmt" ]
then
   echo HERE
   tmin=`head -1 resid.out | awk '{print $1}'`
   tmax=`tail -1 resid.out | awk '{print $1}'`
   dmin=`awk '{print $2}' resid.out | sort -g | head -1 | awk '{print 5*int($1/5)-5}'`
   dmax=`awk '{print $2}' resid.out | sort -g | tail -1 | awk '{print 5*int($1/5)+5}'`
#   R=`awk '{printf "%.3f %.2f\n", $1 +($2-1)/365.25, $3}' data.in | minmax -I1/5`
#   echo $R
   title=`echo $data "detrended, cleaned"`
   awk '{print $1, $2}' resid.out | \
   psxy -JX -R -Ba5Yf1y/a10f5:mm::."$title":WSen -Sc0.03 -X0 -Y-4.0  -O -K >> $psf.ps
fi

#  Save detrended time-series (and cleaned)
cp data.cl data.detrend.cl

######
##  Rebuild adj.in file for future use
#####
rm -f adj.in

Rate=`grep "Rate in units per year" est1.out | tail -1 | awk '{print $6}'`
cat > adj.in <<EOF
R
$Rate
EOF


nsea=`head -1 sea.in`
n=1
while [ "$n" -le "$nsea" ]
do
   ln=`expr $n + 1`
   per=`sed -n ''$ln'p' sea.in`
#   echo $per
   grep "Period of" est1.out | tail -"$nsea" | sed -n ''$n'p' > junk
   camp=`awk '{print $7}' junk`
   samp=`awk '{print $12}' junk`
#   echo $camp $samp
   cat >> adj.in <<EOF
s
$tstart
$per $camp $samp
EOF
   n=$ln
done


noff=0
noff=`head -1 off.in`
noff=`wc -l off.in | awk '{print $1-1}'`
n=1
while [ "$n" -le "$noff" ]
do
   ln=`expr $n + 1`
   toff=`sed -n ''$ln'p' off.in`
   if [ "$netd" = "otr" ]
   then
     off=`grep "Offset number" est1.out | tail -"$noff" | sed -n ''$n'p' | awk '{print $8}'`
   fi
   if [ "$netd" = "gmt" ]
   then
     off=`grep "Offset number" est1.out | tail -"$noff" | sed -n ''$n'p' | awk '{print $7}'`
   fi
   echo OFFSET $toff $off
   cat >> adj.in <<EOF
o
$toff
$off
EOF
   n=$ln
done


nexp=`head -1 exp.in| awk '{print $1}'`

n=1
while [ "$n" -le "$nexp" ]
do
  grep "Exponential number" est1.out | tail -"$nexp" |  sed -n ''$n'p'
  ln=`echo $n | awk '{print 1 + $1*3}'`
  tp=`sed -n ''$ln'p' exp.in`
  ln=`echo $n | awk '{print 1 + $1*3-2}'`
  texp=`sed -n ''$ln'p' exp.in`
  if [ "$nett" = "otr" ]
  then
    amp=`grep "Exponential number" est1.out | tail -"$nexp" |  sed -n ''$n'p' | awk '{print $8}'`
    tau=`grep "Exponential number" est1.out | tail -"$nexp" |  sed -n ''$n'p' | awk '{print $14}'`
  fi
  cat >> adj.in <<EOF
$tp
$texp
$tau $amp
EOF
  n=`expr $n + 1 `
done

nrate=`head -1 rate.in | awk '{print $1}'`

n=1
while [ "$n" -le "$nrate" ]
do
   rate=`grep "Rate change number" est1.out | tail -"$nrate" | sed -n ''$n'p' | awk '{print $12}'`
   ln=`expr $n + 1 `
   trate=`sed -n ''$ln'p' rate.in`
   cat >> adj.in <<EOF
r
$trate
$rate
EOF
   n=`expr $n + 1`
done
#cat adj.in



if [ "$nett" = "otr"  ]
then
  mv $psf.ps "$here"/"$data".ps
#  echo plot in "$here"/"$data".ps
fi
if [ "$nett" = "gmt"  ]
then
  mv $psf.ps "$here"/"$data".ps
# echo plot in "$here"/"$data".ps
fi
mv est_"$data".out "$here"/est_"$data".out
echo "Results can be found in"
ls -1 "$here"/"$data".ps
ls -1 "$here"/est_"$data".out

exit

