#!/bin/bash

#####
#	PARSE COMMAND LINE ARGUMENTS
#####
if [ $# -lt 3 ]
then
    echo "$0: no input files" 1>&2
    echo "Usage: $0 TEMP_FILE DEP_FILE HF_FILE" 1>&2
    echo "  Plot TQTec outputs"
    exit 1
fi

TEMPFILE="$1"
DEPFILE="$2"
HFFILE="$3"
if [ ! -f $TEMPFILE ]
then
    echo "$0: could not find temperature file \"$TEMPFILE\"" 1>&2
    exit 1
fi
if [ ! -f $DEPFILE ]
then
    echo "$0: could not find depth file \"$DEPFILE\"" 1>&2
    exit 1
fi
if [ ! -f $HFFILE ]
then
    echo "$0: could not find heat flow file \"$HFFILE\"" 1>&2
    exit 1
fi

GEOTHERM_FILE=
shift
shift
while [ "$1" != "" ]
do
    case $1 in
        -geotherm) shift; GEOTHERM_FILE=$1;;
    esac
    shift
done


#####
#	SET BOUNDS
#####
# Temperature
TMIN="0"
TMAX=`awk '{if (NR>1) print $0}' $TEMPFILE | gmt gmtinfo -C |\
        awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
        awk '{print $2-$2%20+20}'`

# Depth
ZMIN=`awk '{if (NR>1) print $0}' $DEPFILE | gmt gmtinfo -C |\
      awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $1-$1%10-10}'`
ZMAX=`awk '{if (NR>1) print $0}' $DEPFILE | gmt gmtinfo -C |\
      awk '{for(i=1;i<=NF;i++){printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $2-$2%10+10}'`
ZMAX0="0"

# Time
tMIN=`head -1 $DEPFILE | awk '{print $1}'`
tMAX="0"
dt=`head -1 $DEPFILE | awk '{print $2}'`

# Heat flow
HFMIN=`awk '{if (NR>1) print $0}' $HFFILE | gmt gmtinfo -C |\
      awk '{for (i=1;i<=NF;i++) {printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $1-$1%5-5}'`
HFMAX=`awk '{if (NR>1) print $0}' $HFFILE | gmt gmtinfo -C |\
      awk '{for (i=1;i<=NF;i++) {printf("%f\n"),$i}}' | gmt gmtinfo -C |\
      awk '{print $2-$2%5+5}'`


#####
#	PLOT RESULTS
#####
PSFILE="tqtec_results.ps"

gmt set PS_MEDIA 17ix17i
gmt set MAP_GRID_PEN 0.5p,225,4_4:0

# Colors for horizons
NCOL=$(sed -ne "2p" $DEPFILE | awk '{print NF}')
NCOL2=$(echo $NCOL 1.1 | awk '{print $1*$2}')
gmt makecpt -Cmagma -T0/$NCOL2/1 | awk '{if(NF==5){print $2}}' > color_list.tmp


# Depth vs. time
echo "$0: plotting depth vs. time"

LIMS="-R$tMIN/$tMAX/$ZMIN/$ZMAX"
PROJ="-JX4i/3i -P"

echo 0 0 | gmt psxy $PROJ $LIMS -K -X1.5i -Y1.5i > $PSFILE

gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya10g5+l"Depth (km)" -BWeSn -K -O >> $PSFILE

for COL in $(seq 1 $NCOL)
do
    COLOR=$(sed -ne "${COL}p" color_list.tmp)
    if [ $COL -eq 1 ]
    then
        awk 'BEGIN{g=1}{
            if (NR == 1) {print "> -W3p,'$COLOR'"}
            if (NR > 1) {
                if (g==1 && $'$COL'>0) {
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
                if (g==1 && $'$COL'>0) {
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

echo 0.05 0.05 10,2 LB $(date) | gmt pstext -JX1i -R0/1/0/1 -F+f+j -N -K -O >> $PSFILE



# Temperature vs. time
echo "$0: plotting temperature vs. time"

LIMS="-R$tMIN/$tMAX/$TMIN/$TMAX"

echo 0 0 | gmt psxy $PROJ $LIMS -Y4i -K -O >> $PSFILE

gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya20g20+l"Temperature (\260C)" -BWeSn -K -O >> $PSFILE

for COL in $(seq 1 $NCOL)
do
    COLOR=$(sed -ne "${COL}p" color_list.tmp)
    if [ $COL -eq 1 ]
    then
        awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$'$COL'}' $TEMPFILE |\
          gmt psxy $PROJ $LIMS -W3p,$COLOR -K -O >> $PSFILE
    else
        awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$'$COL'}' $TEMPFILE |\
          gmt psxy $PROJ $LIMS -W1p,$COLOR -K -O >> $PSFILE
    fi
done



# Surface heat flow vs. time
echo "$0: plotting surface heat flow vs. time"

LIMS="-R$tMIN/$tMAX/$HFMIN/$HFMAX"

echo 0 0 | gmt psxy $PROJ $LIMS -Y4i -K -O >> $PSFILE

gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya5g5+l"Surface Heat Flow (W/m*K)" -BWeSn -K -O >> $PSFILE

awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$1}' $HFFILE |\
    gmt psxy $PROJ $LIMS -W1p -K -O >> $PSFILE



# Depth vs. temperature
echo "$0: plotting depth vs. temperature"
LIMS="-R$TMIN/$TMAX/$ZMIN/$ZMAX0"
PROJ="-JX5i/5i"

echo 0 0 | gmt psxy $PROJ $LIMS -X5.5i -Y-4i -K -O >> $PSFILE

gmt psbasemap $PROJ $LIMS -Bxa20g20+l"Temperature (C)" -Bya10g5+l"Depth (km)" -BWeSn -K -O >> $PSFILE

if [ "$GEOTHERM_FILE" != "" -a "$PLOT_GEOTHERM_DEP_VS_TEMP" == "Y" ]
then
    awk '{
        if (/>/) {
            if ($3==0) {
                print "> -W1p"
            } else {
                print "> -W0.5p,105"
            }
        } else {
            print $1,-$2
        }
    }' $GEOTHERM_FILE |\
        gmt psxy $PROJ $LIMS -K -O >> $PSFILE
fi

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


# Temperature-depth contours over time
if [ "$GEOTHERM_FILE" != "" ]
then
    echo "$0: plotting temperature contours over time"
    LIMS="-R$tMIN/$tMAX/$ZMIN/$ZMAX"
    PROJ="-JX4i/3i -P"

    echo 0 0 | gmt psxy $PROJ $LIMS -K -O -Y-4i >> $PSFILE

    gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya10g5+l"Depth (km)" -BWeSn -K -O >> $PSFILE

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
fi


echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

gmt psconvert $PSFILE -Tg -A
echo "Created file $(basename $PSFILE .ps).png"

rm $PSFILE
rm gmt.*
rm color_list.tmp
