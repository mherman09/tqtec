#!/bin/bash

if [ $# -ne 1 ]
then
    echo "$0: no input file" 1>&2
    echo "Usage: $0 TQTec_depth_file" 1>&2
    echo "  Plot the depths of units tracked in TQTec" 1>&2
    echo "  Eroded units plotted at ground level" 1>&2
    exit 1
fi

DEPFILE="$1"
if [ ! -f $DEPFILE ]
then
    echo "$0: could not find file \"$DEPFILE\"" 1>&2
    exit 1
fi


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
gmt set PS_MEDIA 11ix11i
gmt set MAP_GRID_PEN 0.5p,225,4_4:0

LIMS="-R$tMIN/$tMAX/$ZMIN/$ZMAX"
PROJ="-JX8i/6i"
PSFILE="dep0_vs_time.ps"


#####
# Plot depth vs. time
#####
gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya10g5+l"Depth (km)" -BWeSn -X1.35i -Y1.25i -K > $PSFILE

NCOL=$(sed -ne "2p" $DEPFILE | awk '{print NF}')
NCOL2=$(echo $NCOL 1.1 | awk '{print $1*$2}')
gmt makecpt -Cmagma -T0/$NCOL2/1 | awk '{if(NF==5){print $2}}' > color_list.tmp
for COL in $(seq 1 $NCOL)
do
    COLOR=$(sed -ne "${COL}p" color_list.tmp)
    if [ $COL -eq 1 ]
    then
        awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$'$COL'}' $DEPFILE |\
            awk '{if ($2 > 0) print $1,0;else print $0}' |\
            gmt psxy $PROJ $LIMS -W3p,$COLOR -K -O >> $PSFILE
    else
        awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$'$COL'}' $DEPFILE |\
            awk '{if ($2 > 0) print $1,0;else print $0}' |\
            gmt psxy $PROJ $LIMS -W1p,$COLOR -K -O >> $PSFILE
    fi
done

gmt psxy $PROJ $LIMS -W0.5p,105,12_4:0 -K -O >> $PSFILE << EOF
 $tMIN 0
 $tMAX 0
EOF

echo 0.05 0.05 10,2 LB $(date) | gmt pstext -JX1i -R0/1/0/1 -F+f+j -N -K -O >> $PSFILE

echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

gmt psconvert $PSFILE -Tg -A
echo "Created file $(basename $PSFILE .ps).png"

rm $PSFILE
rm gmt.*
rm color_list.tmp
