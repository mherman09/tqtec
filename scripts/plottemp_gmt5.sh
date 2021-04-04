#!/bin/sh

if [ $# -ne 1 ]
then
    echo "Error: No input file!"
    echo "Usage: plottemp.sh [TQTec_temp_file]"
    echo "  Plot the temperatures of units tracked in TQTec."
    exit
fi

TEMPFILE="$1"

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
LIMS="-R$TIMEMIN/$TIMEMAX/$TEMPMIN/$TEMPMAX"
PROJ="-JX6i"
PSFILE="Tvst.ps"

#####
# Plot temperature vs. time
#####
gmt psbasemap $PROJ $LIMS -Bxa10+l"Time (Ma)" -Bya20+l"Temperature (\260C)" -BWeSn -K > $PSFILE

for COL in 1 2 3 4 5 6 7 8 9 10
do
    if [ $COL -eq 1 ]
    then
        awk '{if (NR > 1) print '$TIMEMIN'+(NR-1)*'$dt',$'$COL'}' $TEMPFILE |\
          gmt psxy $PROJ $LIMS -W3p -K -O >> $PSFILE
    else
        awk '{if (NR > 1) print '$TIMEMIN'+(NR-1)*'$dt',$'$COL'}' $TEMPFILE |\
          gmt psxy $PROJ $LIMS -W1p -K -O >> $PSFILE
    fi
done

echo "Created file $PSFILE"
ps2pdf $PSFILE
