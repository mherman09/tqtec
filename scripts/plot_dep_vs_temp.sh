#!/bin/sh

if [ $# -lt 2 ]
then
    echo "$0: no input file" 1>&2
    echo "Usage: $0 TQTec_temp_file TQTec_depth_file [-init temp_init_file]" 1>&2
    echo "  Plot the temperature-depth paths of tracked units over time"
    exit 1
fi

TEMPFILE="$1"
DEPFILE="$2"
if [ ! -f $TEMPFILE ]
then
    echo "$0: could not find file \"$TEMPFILE\"" 1>&2
    exit 1
fi
if [ ! -f $DEPFILE ]
then
    echo "$0: could not find file \"$DEPFILE\"" 1>&2
    exit 1
fi

TEMP_INIT_FILE=
shift
shift
while [ "$1" != "" ]
do
    case $1 in
        -init) shift; TEMP_INIT_FILE=$1;;
    esac
    shift
done

#####
# Temperature bounds
#####
TMIN="0"
TMAX=`awk '{if (NR>1) print $0}' $TEMPFILE | gmt gmtinfo -C |\
        awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
        awk '{print $2-$2%20+20}'`

#####
# Depth bounds
#####
ZMIN=`awk '{if (NR>1) print $0}' $DEPFILE | gmt gmtinfo -C |\
        awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
        awk '{print $1-$1%10-10}'`
ZMAX="0"

#####
# GMT variables
#####
LIMS="-R$TMIN/$TMAX/$ZMIN/$ZMAX"
PROJ="-JX6i/6i"
PSFILE="dep_vs_temp.ps"

#####
# Plot depth-temperature
#####
gmt psbasemap $PROJ $LIMS -Bxa20+l"Temperature (C)" -Bya5+l"Depth (km)" -BWeSn -K > $PSFILE

if [ "$TEMP_INIT_FILE" != "" ]
then
    awk '{print $1,-$2}' $TEMP_INIT_FILE |\
        gmt psxy $PROJ $LIMS -W1p -K -O >> $PSFILE
fi

NCOL=$(sed -ne "2p" $TEMPFILE | awk '{print NF}')
NCOL2=$(echo $NCOL 1.1 | awk '{print $1*$2}')
gmt makecpt -Cmagma -T0/$NCOL2/1 | awk '{if(NF==5){print $2}}' > color_list.tmp
for i in $(seq 1 $NCOL)
do
    COLOR=$(sed -ne "${i}p" color_list.tmp | awk '{print $1}')
    # Plot curves
    paste $TEMPFILE $DEPFILE |\
        awk '{if (NR > 1) print $'$i',$('$i'+10)}' |\
        gmt psxy $PROJ $LIMS -W1p,$COLOR -K -O >> $PSFILE
    # Large dot at 0 Ma
    paste $TEMPFILE $DEPFILE |\
        awk '{if (NR == 2) print $'$i',$('$i'+10)}' |\
        gmt psxy $PROJ $LIMS -Sc0.10i -G$COLOR -W0.5p -K -O >> $PSFILE
    # Dots every 10 Ma
    paste $TEMPFILE $DEPFILE |\
        awk '{if (NR>1 && ((NR-1)*0.01)%10 == 0) print $'$i',$('$i'+10)}' |\
        gmt psxy $PROJ $LIMS -Sc0.05i -G$COLOR -W0.5p -K -O >> $PSFILE
done

echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

gmt psconvert $PSFILE -Tg -A
echo "Created file $(basename $PSFILE .ps).png"
