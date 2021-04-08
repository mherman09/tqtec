#!/bin/bash

if [ $# -ne 1 ]
then
    echo "$0: no input file" 1>&2
    echo "Usage: $0 TQTec_depth_file" 1>&2
    echo "  Plot the depths of units tracked in TQTec" 1>&2
    echo "  Eroded units plotted above ground level" 1>&2
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
PSFILE="dep_vs_time.ps"


#####
# Plot depth vs. time
#####
gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya10g5+l"Depth (km)" -BWESn -X1.35i -Y1.25i -K > $PSFILE

NCOL=$(sed -ne "2p" $DEPFILE | awk '{print NF}')
NCOL2=$(echo $NCOL 1.1 | awk '{print $1*$2}')
gmt makecpt -Cmagma -T0/$NCOL2/1 | awk '{if(NF==5){print $2}}' > color_list.tmp
for COL in $(seq 1 $NCOL)
do
    COLOR=$(sed -ne "${COL}p" color_list.tmp)
    if [ $COL -eq 1 ]
    then
        awk 'BEGIN{g=1}{
            if (NR == 1) {print "> -W3p,'$COLOR'"}
            if (NR > 1) {
                if (g==1 && $'$COL'>=0) {
                    g = 0
                    print "> -W3p,'$COLOR'@80"
                }
                if (g==0 && $'$COL'<=0) {
                    g = 1
                    print "> -W3p,'$COLOR'"
                }
                print '$tMIN'+(NR-1)*'$dt',$'$COL'
            }
        }' $DEPFILE |\
            gmt psxy $PROJ $LIMS -K -O >> $PSFILE
    else
        awk 'BEGIN{g=1}{
            if (NR == 1) {print "> -W1p,'$COLOR'"}
            if (NR > 1) {
                if (g==1 && $'$COL'>=0) {
                    g = 0
                    print "> -W1p,'$COLOR'@80"
                }
                if (g==0 && $'$COL'<=0) {
                    g = 1
                    print "> -W1p,'$COLOR'"
                }
                print '$tMIN'+(NR-1)*'$dt',$'$COL'
            }
        }' $DEPFILE |\
            gmt psxy $PROJ $LIMS -K -O >> $PSFILE
    fi
done

gmt psxy $PROJ $LIMS -W1p,105,12_4:0 -K -O >> $PSFILE << EOF
 $tMIN 0
 $tMAX 0
EOF

echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

gmt psconvert $PSFILE -Tg -A
echo "Created file $(basename $PSFILE .ps).png"
