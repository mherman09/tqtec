#!/bin/bash

if [ $# -ne 1 ]
then
    echo "$0: no input file" 1>&2
    echo "Usage: $0 TQTec_temp_file" 1>&2
    echo "  Plot the temperatures of units tracked in TQTec" 1>&2
    exit 1
fi

TEMPFILE="$1"
if [ ! -f $TEMPFILE ]
then
    echo "$0: could not find file \"$TEMPFILE\"" 1>&2
    exit 1
fi

#####
# Temperature bounds
#####
TEMPMIN="0"
TEMPMAX=`awk '{if (NR>1) print $0}' $TEMPFILE | gmt gmtinfo -C |\
      awk '{for (i=1;i<=NF;i++) {printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $2-$2%20+20}'`

#####
# Time bounds
#####
TIMEMIN=`head -1 $TEMPFILE | awk '{print $1}'`
TIMEMAX="0"
dt=`head -1 $TEMPFILE | awk '{print $2}'`

#####
# GMT variables
#####
gmt set PS_MEDIA 11ix11i

LIMS="-R$TIMEMIN/$TIMEMAX/$TEMPMIN/$TEMPMAX"
PROJ="-JX6i/6i"
PSFILE="temp_vs_time.ps"

#####
# Plot temperature vs. time
#####
gmt psbasemap $PROJ $LIMS -Bxa10+l"Time (Ma)" -Bya20+l"Temperature (\260C)" -BWeSn -K > $PSFILE

NCOL=$(sed -ne "2p" $TEMPFILE | awk '{print NF}')
NCOL2=$(echo $NCOL 1.1 | awk '{print $1*$2}')
gmt makecpt -Cmagma -T0/$NCOL2/1 | awk '{if(NF==5){print $2}}' > color_list.tmp
for COL in $(seq 1 $NCOL)
do
    COLOR=$(sed -ne "${COL}p" color_list.tmp)
    if [ $COL -eq 1 ]
    then
        awk '{if (NR > 1) print '$TIMEMIN'+(NR-1)*'$dt',$'$COL'}' $TEMPFILE |\
          gmt psxy $PROJ $LIMS -W3p,$COLOR -K -O >> $PSFILE
    else
        awk '{if (NR > 1) print '$TIMEMIN'+(NR-1)*'$dt',$'$COL'}' $TEMPFILE |\
          gmt psxy $PROJ $LIMS -W1p,$COLOR -K -O >> $PSFILE
    fi
done

echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

gmt psconvert $PSFILE -Tg -A
echo "Created file $(basename $PSFILE .ps).png"
