#!/bin/bash

#####
#	PARSE COMMAND LINE
#####
function usage () {
    echo "Usage: $0 HF_FILE [-timing ACTION_TIME_FILE] [-n NPTAVG] [-o BASENAME]" 1>&2
    echo "" 1>&2
    echo "HF_FILE                     Surface heat flow versus time file generated by readtqtec" 1>&2
    echo "-timing ACTION_TIME_FILE    Timing of tectonic actions generated by tqtec -timing" 1>&2
    echo "-n NPTAVG                   Superimpose n-point average of heat flow" 1>&2
    echo "-o BASENAME                 Basename for output file (default: hf_vs_time)" 1>&2
    exit 1
}

# Required heat flow vs. time file
HFFILE="$1"
if [ ! -f $HFFILE ]
then
    echo "$0: could not find heat flow file \"$HFFILE\"" 1>&2
    exit 1
fi
shift

# Optional arguments
TIMING_FILE=
NPTAVG=
BASENAME=hf_vs_time
while [ "$1" != "" ]
do
    case $1 in
        -timing) shift; TIMING_FILE=$1;;
        -n) shift; NPTAVG=$1;;
        -o) shift; BASENAME=$1;;
    esac
    shift
done


#####
#	SET BOUNDS
#####
# Heat flow
HFMIN=`awk '{if (NR>1) print $0}' $HFFILE | gmt gmtinfo -C |\
      awk '{for (i=1;i<=NF;i++) {printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $1-$1%5-5}'`
HFMAX=`awk '{if (NR>1) print $0}' $HFFILE | gmt gmtinfo -C |\
      awk '{for (i=1;i<=NF;i++) {printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $2-$2%5+5}'`

# Time
tMIN=`head -1 $HFFILE | awk '{print $1}'`
tMAX="0"
dt=`head -1 $HFFILE | awk '{print $2}'`

# Tectonic action timing
function plot_tectonic_timing () {
    YMIN=$(echo $LIMS | sed -e "s/-R//" | awk -F/ '{print $3}')
    YMAX=$(echo $LIMS | sed -e "s/-R//" | awk -F/ '{print $4}')
    awk 'BEGIN{print ">"}{
        if ($1=="burial" || $1=="uplift") {
            print ">"
            print '$tMIN'+$3,'$YMIN'
            print '$tMIN'+$4,'$YMIN'
            print '$tMIN'+$4,'$YMAX'
            print '$tMIN'+$3,'$YMAX'
            print '$tMIN'+$3,'$YMIN'
        }
    }' $TIMING_FILE |\
        gmt psxy $PROJ $LIMS -G245 -K -O >> $PSFILE
    awk 'BEGIN{print ">"}{
        if ($1=="thrust") {
            print ">"
            print '$tMIN'+$3,'$YMIN'
            print '$tMIN'+$3,'$YMAX'
        }
    }' $TIMING_FILE |\
        gmt psxy $PROJ $LIMS -W3p,225 -K -O >> $PSFILE
}


#####
# GMT variables
#####

gmt set PS_MEDIA 17ix17i
gmt set MAP_GRID_PEN 0.5p,225,4_4:0

LIMS="-R$tMIN/$tMAX/$HFMIN/$HFMAX"
tMIN_MODEL=`echo $tMIN | awk '{print -$1}'`
tMAX_MODEL=`echo $tMAX | awk '{print -$1}'`
LIMS_MODEL="-R$tMAX_MODEL/$tMIN_MODEL/$HFMIN/$HFMAX"

AGE_TIKS=`echo $LIMS | sed -e "s/-R//" | awk -F"/" '{print $1,$2}' |\
          awk '{
              time = $2-$1
              if (time>=50) {dt=10}
              else if (time>=20) {dt=5}
              else if (time>=10) {dt=2}
              else {dt=1}
              print dt
          }'`
HF_TIKS=`echo $LIMS | sed -e "s/-R//" | awk -F"/" '{print $3,$4}' |\
         awk '{
            hf = $2-$1
            if (hf>=200) {dhf=50}
            else if (hf>=100) {dhf=20}
            else if (hf>=50) {dhf=10}
            else if (hf>=20) {dhf=5}
            else if (hf>=10) {dhf=2}
            else if (hf>=5) {dhf=1}
            else {dhf=0.5}
            print dhf
          }'`

PROJ="-JX8i/6i"

PSFILE="${BASENAME}.ps"


# Initialize plot
gmt psxy -T -K > $PSFILE

if [ "$TIMING_FILE" != "" ]; then plot_tectonic_timing; fi

gmt psbasemap $PROJ $LIMS -Bxa${AGE_TIKS}g${AGE_TIKS}+l"Time Before Present (Ma)" -Bya${HF_TIKS}g${HF_TIKS}+l"Heat Flow (mW/m@+2@+)" -BWeS -K -O >> $PSFILE
gmt psbasemap $PROJ $LIMS_MODEL -Bxa${AGE_TIKS}+l"Model Time (Ma)" -BN -K -O >> $PSFILE

# Surface heat flow versus time
awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$1}' $HFFILE |\
    gmt psxy $PROJ $LIMS -W1p,black -K -O >> $PSFILE

# N-point average heat flow versus time
if [ "$NPTAVG" != "" ]
then
    NPTS=`wc $HFFILE | awk '{print $1-1}'`
    awk '{
        if (NR>1) {
            hf[NR-1] = $1
        }
    } END {
        sum = 0
        for (i=1;i<NR;i++) {
            if (i<int('$NPTAVG'/2)) {
                print hf[i]
            } else if (i>'$NPTS'-int('$NPTAVG'/2)) {
                print hf[i]
            } else {
                sum = 0
                for (j=-int('$NPTAVG'/2);j<=int('$NPTAVG'/2)-1;j++) {
                    sum = sum + hf[i+j]/'$NPTAVG'
                }
                print sum
            }
        }
    }' $HFFILE |
        awk '{print '$tMIN'+(NR-1)*'$dt',$1}' |\
        gmt psxy $PROJ $LIMS -W2p,155/55/25 -K -O >> $PSFILE
fi

# Date and time
echo 12,0 LB $(date) | gmt pstext $PROJ $LIMS -F+f+j+cBL -D0.05i -K -O >> $PSFILE

# Finalize
gmt psxy -T -O >> $PSFILE

gmt psconvert $PSFILE -Tf -A
echo "Created file $(basename $PSFILE .ps).pdf"

rm $PSFILE
rm gmt.*
