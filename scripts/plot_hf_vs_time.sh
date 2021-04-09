#!/bin/bash

if [ $# -ne 1 ]
then
    echo "$0: no input file" 1>&2
    echo "Usage: $0 TQTec_hf_file" 1>&2
    echo "  Plot the surface heat flow over time from TQTec model" 1>&2
    exit 1
fi

HFFILE="$1"
if [ ! -f $HFFILE ]
then
    echo "$0: could not find file \"$HFFILE\"" 1>&2
    exit 1
fi

#####
# Heat flow bounds
#####
HFMIN=`awk '{if (NR>1) print $0}' $HFFILE | gmt gmtinfo -C |\
      awk '{for (i=1;i<=NF;i++) {printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $1-$1%5-5}'`
HFMAX=`awk '{if (NR>1) print $0}' $HFFILE | gmt gmtinfo -C |\
      awk '{for (i=1;i<=NF;i++) {printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $2-$2%5+5}'`

#####
# Time bounds
#####
TIMEMIN=`head -1 $HFFILE | awk '{print $1}'`
TIMEMAX="0"
dt=`head -1 $HFFILE | awk '{print $2}'`

#####
# GMT variables
#####
gmt set PS_MEDIA 11ix11i
gmt set MAP_GRID_PEN 0.5p,225,4_4:0

LIMS="-R$TIMEMIN/$TIMEMAX/$HFMIN/$HFMAX"
PROJ="-JX6i/6i"
PSFILE="hf_vs_time.ps"

#####
# Plot heat flow vs. time
#####
gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya5g5+l"Surface Heat Flow (W/m*K)" -BWeSn -K > $PSFILE

awk '{if (NR > 1) print '$TIMEMIN'+(NR-1)*'$dt',$1}' $HFFILE |\
    gmt psxy $PROJ $LIMS -W1p -K -O >> $PSFILE

echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

gmt psconvert $PSFILE -Tg -A
echo "Created file $(basename $PSFILE .ps).png"

rm $PSFILE
rm gmt.*
