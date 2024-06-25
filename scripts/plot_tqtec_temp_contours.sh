#!/bin/bash

#####
#	PARSE COMMAND LINE
#####
function usage () {
    echo "Usage: $0 DEP_FILE GEOTHERM_FILE [-timing ACTION_TIME_FILE] [-o BASENAME]" 1>&2
    echo "" 1>&2
    echo "DEP_FILE                    Depth versus time file generated by readtqtec" 1>&2
    echo "GEOTHERM_FILE               Geotherm generated by tqtec -geotherm" 1>&2
    echo "-timing ACTION_TIME_FILE    Timing of tectonic actions generated by tqtec -timing" 1>&2
    echo "-o BASENAME                 Basename for output file (default: dep_vs_time)" 1>&2
    exit 1
}

# Required depth vs. time file
DEPFILE="$1"
GEOTHERM_FILE="$2"
if [ ! -f $DEPFILE ]
then
    echo "$0: could not find depth file \"$DEPFILE\"" 1>&2
    exit 1
fi
if [ ! -f $GEOTHERM_FILE ]
then
    echo "$0: could not find depth file \"$GEOTHERM_FILE\"" 1>&2
    exit 1
fi
shift

# Optional arguments
TIMING_FILE=
BASENAME=temp_contours
while [ "$1" != "" ]
do
    case $1 in
        -timing) shift; TIMING_FILE=$1;;
        -o) shift; BASENAME=$1;;
    esac
    shift
done


#####
#	SET BOUNDS
#####
# Depth
ZMIN=`awk '{if (NR>1) print $0}' $DEPFILE | gmt gmtinfo -C |\
      awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $1-$1%10-10}'`
ZMAX=`awk '{if (NR>1) print $0}' $DEPFILE | gmt gmtinfo -C |\
      awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $2-$2%10+10}'`

# Time
tMIN=`head -1 $DEPFILE | awk '{print $1}'`
tMAX="0"
dt=`head -1 $DEPFILE | awk '{print $2}'`

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
#	PLOT TEMPERATURE CONTOURS
#####

# GMT variables
gmt set PS_MEDIA 17ix17i
gmt set MAP_GRID_PEN 0.5p,225,4_4:0

LIMS="-R$tMIN/$tMAX/$ZMIN/$ZMAX"
tMIN_MODEL=`echo $tMIN | awk '{print -$1}'`
tMAX_MODEL=`echo $tMAX | awk '{print -$1}'`
LIMS_MODEL="-R$tMAX_MODEL/$tMIN_MODEL/$ZMIN/$ZMAX"

AGE_TIKS=`echo $LIMS | sed -e "s/-R//" | awk -F"/" '{print $1,$2}' |\
          awk '{
              time = $2-$1
              if (time>=50) {dt=10}
              else if (time>=20) {dt=5}
              else if (time>=10) {dt=2}
              else {dt=1}
              print dt
          }'`
DEP_TIKS=`echo $LIMS | sed -e "s/-R//" | awk -F"/" '{print $3,$4}' |\
          awk '{
              relief = $2-$1
              if (relief>=50) {dz=10}
              else if (relief>=20) {dz=5}
              else if (relief>=10) {dz=2}
              else if (relief>=5) {dz=1}
              else {dz=0.5}
              print dz
          }'`

PROJ="-JX8i/6i"

PSFILE="${BASENAME}.ps"


# Colors for each horizon
NCOL=$(sed -ne "2p" $DEPFILE | awk '{print NF}')
NCOL2=$(echo $NCOL 1.1 | awk '{print $1*$2}')
gmt makecpt -Cmagma -T0/$NCOL2/1 | awk '{if(NF==5){print $2}}' > color_list.tmp


# Initialize plot
gmt psxy -T -K > $PSFILE

# Tectonic timing
if [ "$TIMING_FILE" != "" ]; then plot_tectonic_timing; fi

# Basemap
gmt psbasemap $PROJ $LIMS -Bxa${AGE_TIKS}g${AGE_TIKS}+l"Time Before Present (Ma)" -Bya${DEP_TIKS}g${DEP_TIKS}+l"Depth (km)" -BWeS -K -O >> $PSFILE
gmt psbasemap $PROJ $LIMS_MODEL -Bxa${AGE_TIKS}+l"Model Time (Ma)" -BN -K -O >> $PSFILE

# Depth versus time
DEP_DASH="2_1:0"
for COL in $(seq 1 $NCOL)
do
    COLOR=$(sed -ne "${COL}p" color_list.tmp)
    if [ $COL -eq 1 ]
    then
        awk 'BEGIN{g=1}{
            if (NR == 1) {print "> -W1.5p,'$COLOR','$DEP_DASH'"}
            if (NR > 1) {
                if (g==1 && $'$COL'>0) {
                    g = 0
                    print "> -W1.5p,'$COLOR'@80,'$DEP_DASH'"
                }
                if (g==0 && $'$COL'<=0) {
                    g = 1
                    print "> -W1.5p,'$COLOR','$DEP_DASH'"
                }
                print '$tMIN'+(NR-1)*'$dt',$'$COL'
            }
        }' $DEPFILE |\
            gmt psxy $PROJ $LIMS -K -O >> $PSFILE
    else
        awk 'BEGIN{g=1}{
            if (NR == 1) {print "> -W1p,'$COLOR','$DEP_DASH'"}
            if (NR > 1) {
                if (g==1 && $'$COL'>0) {
                    g = 0
                    print "> -W1p,'$COLOR'@80,'$DEP_DASH'"
                }
                if (g==0 && $'$COL'<=0) {
                    g = 1
                    print "> -W1p,'$COLOR','$DEP_DASH'"
                }
                print '$tMIN'+(NR-1)*'$dt',$'$COL'
            }
        }' $DEPFILE |\
            gmt psxy $PROJ $LIMS -K -O >> $PSFILE
    fi
done

# Temperature contours
t1=$(echo $tMIN $tMAX | awk '{print $1+($2-$1)/20}')
t2=$(echo $tMIN $tMAX | awk '{print $1+($2-$1)/20}')
t3=$(echo $tMIN $tMAX | awk '{print $2-($2-$1)/20}')
t4=$(echo $tMIN $tMAX | awk '{print $2-($2-$1)/20}')
Z1=$ZMIN
Z2=$ZMAX
awk '{
    if (/>/) {
        time = '$tMIN' + $3*'$dt'/2
    } else {
        print time,-$2,$1
    }
}' $GEOTHERM_FILE |\
    gmt pscontour $PROJ $LIMS -Wa1p -Wc0.5p -C10 -A50+f7p+u"\260C" -Gl$t1/$Z1/$t2/$Z2,$t3/$Z1/$t4/$Z2 -K -O >> $PSFILE

echo $tMIN $tMAX | awk '{print $1,0;print $2,0}' | gmt psxy $PROJ $LIMS -W1p,105,12_4:0 -K -O >> $PSFILE

# Date and time
echo 12,0 LB $(date) | gmt pstext $PROJ $LIMS -F+f+j+cBL -D0.05i -K -O >> $PSFILE

# Finalize
gmt psxy -T -O >> $PSFILE

gmt psconvert $PSFILE -Tf -A
echo "Created file $(basename $PSFILE .ps).pdf"

rm $PSFILE
rm gmt.*
rm color_list.tmp