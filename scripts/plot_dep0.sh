#!/bin/sh

if [ $# -ne 1 ]
then
    echo "Error: No input file!"
    echo "Usage: plotdep0.sh [TQTec_depth_file]"
    echo "  Plot the depths of units tracked in TQTec."
    echo "  Eroded units plotted at ground level."
    exit
fi

DEPFILE="$1"

#####
# Depth bounds
#####
ZMIN=`awk '{if (NR>1) print $0}' $DEPFILE | gmt gmtinfo -C |\
      awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $1-$1%10-10}'`
ZMAX=`awk '{if (NR>1) print $0}' $DEPFILE | gmt gmtinfo -C |\
      awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $2-$2%10+10}'`

#####
# Time bounds
#####
tMIN=`head -1 $DEPFILE | awk '{print $1}'`
tMAX="0"
dt=`head -1 $DEPFILE | awk '{print $2}'`

#####
# GMT variables
#####
LIMS="-R$tMIN/$tMAX/$ZMIN/$ZMAX"
PROJ="-JX8.3i/6.5i"
PSFILE="Z0vst.ps"

gmt set MAP_GRID_PEN 0.5p,155/155/155,4_4:0

#####
# Plot depth vs. time
#####
gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya10g10+l"Depth (km)" -BWESn -X1.35i -Y1.25i -K > $PSFILE

for COL in 1 2 3 4 5 6 7 8 9 10
do
    if [ $COL -eq 1 ]
    then
        awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$'$COL'}' $DEPFILE |\
          awk '{if ($2 > 0) print $1,0;else print $0}' |\
          gmt psxy $PROJ $LIMS -W3p -K -O >> $PSFILE
    else
        awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$'$COL'}' $DEPFILE |\
          awk '{if ($2 > 0) print $1,0;else print $0}' |\
          gmt psxy $PROJ $LIMS -W1p -K -O >> $PSFILE
    fi  
done

gmt psxy $PROJ $LIMS -W1p,0/0/0,4_4:0 -K -O >> $PSFILE << EOF
$tMIN 0
$tMAX 0
EOF

echo "Created file $PSFILE"
ps2pdf $PSFILE
