#!/bin/bash

#####
#	PARSE COMMAND LINE ARGUMENTS
#####
function usage () {
    echo "Usage: $0 TEMP_FILE  DEP_FILE  HF_FILE  [-geotherm GEOTHERM_FILE [-geotherm:time TIME]] [-timing ACTION_TIME_FILE] [-closure CLOSURE_FILE]" 1>&2
    echo "" 1>&2
    echo "TEMP_FILE                   Temperature versus time file generated by readtqtec" 1>&2
    echo "DEP_FILE                    Depth versus time file generated by readtqtec" 1>&2
    echo "HF_FILE                     Surface heat flow versus time file generated by readtqtec" 1>&2
    echo "-geotherm GEOTHERM_FILE     Geotherm generated by tqtec -geotherm (default: plot geotherm at start of model)" 1>&2
    echo "-geotherm:time TIME         Time to plot geotherm in Ma since start of model (option can be repeated)" 1>&2
    echo "-timing ACTION_TIME_FILE    Timing of tectonic actions generated by tqtec -timing" 1>&2
    echo "-closure CLOSURE_FILE       Plot closure on temp-time plot from readtqtec -closure" 1>&2
    echo "-closure:LIMS -RX1/X2/Y1/Y2 Set closure temperature limites" 1>&2
    exit 1
}

if [ $# -lt 3 ]
then
    echo "$0: three input files must be defined" 1>&2
    usage
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
GEOTHERM_PLOT_TIME=0
TIMING_FILE=
CLOSURE_FILE=
CLOSURE_LIMS=
shift
shift
while [ "$1" != "" ]
do
    case $1 in
        -geotherm) shift; GEOTHERM_FILE=$1;;
        -geotherm:time) shift; GEOTHERM_PLOT_TIME="$GEOTHERM_PLOT_TIME $1";;
        -timing) shift; TIMING_FILE=$1;;
        -closure) shift; CLOSURE_FILE=$1;;
        -closure:LIMS) shift; CLOSURE_LIMS=$1;;
    esac
    shift
done
if [ "$GEOTHERM_FILE" != "" ]
then
    if [ ! -f $GEOTHERM_FILE ]
    then
        echo "$0: could not find geotherm file \"$GEOTHERM_FILE\"...not plotting" 1>&2
        GEOTHERM_FILE=
    fi
fi
if [ "$TIMING_FILE" != "" ]
then
    if [ ! -f $TIMING_FILE ]
    then
        echo "$0: could not find tectonic action timing file \"$TIMING_FILE\"...not plotting" 1>&2
        TIMING_FILE=
    fi
fi
if [ "$CLOSURE_FILE" != "" ]
then
    if [ ! -f $CLOSURE_FILE ]
    then
        echo "$0: could not find closure temperature file \"$CLOSURE_FILE\"...not plotting" 1>&2
        CLOSURE_FILE=
    fi
fi


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

# Tectonic action function
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
tMIN_MODEL=`echo $tMIN | awk '{print -$1}'`
tMAX_MODEL=`echo $tMAX | awk '{print -$1}'`
LIMS_MODEL="-R$tMAX_MODEL/$tMIN_MODEL/$ZMIN/$ZMAX"
PROJ="-JX4i/3i -P"

echo 0 0 | gmt psxy $PROJ $LIMS -K -X1.5i -Y1.5i > $PSFILE

if [ "$TIMING_FILE" != "" ]; then plot_tectonic_timing; fi

gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time Before Present (Ma)" -Bya10g5+l"Depth (km)" -BWeS -K -O >> $PSFILE
gmt psbasemap $PROJ $LIMS_MODEL -Bxa10+l"Model Time (Ma)" -BN -K -O >> $PSFILE

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
LIMS_MODEL="-R$tMAX_MODEL/$tMIN_MODEL/$TMIN/$TMAX"

echo 0 0 | gmt psxy $PROJ $LIMS -Y4.5i -K -O >> $PSFILE

if [ "$TIMING_FILE" != "" ]; then plot_tectonic_timing; fi

gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time Before Present (Ma)" -Bya20g20+l"Temperature (\260C)" -BWeS -K -O >> $PSFILE
gmt psbasemap $PROJ $LIMS_MODEL -Bxa10+l"Model Time (Ma)" -BN -K -O >> $PSFILE

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

if [ "$CLOSURE_FILE" != "" ]
then
    awk '{if(/>/){temp=$2}else{print '$tMIN'+$1,temp}}' $CLOSURE_FILE |\
        gmt psxy $PROJ $LIMS -Sc0.05i -W1p -K -O >> $PSFILE
fi



# Surface heat flow vs. time
echo "$0: plotting surface heat flow vs. time"

LIMS="-R$tMIN/$tMAX/$HFMIN/$HFMAX"

echo 0 0 | gmt psxy $PROJ $LIMS -Y4.5i -K -O >> $PSFILE

if [ "$TIMING_FILE" != "" ]; then plot_tectonic_timing; fi

gmt psbasemap $PROJ $LIMS -Bxa10g10+l"Time (Ma)" -Bya5g5+l"Surface Heat Flow (mW/m@+2@+)" -BWeSn -K -O >> $PSFILE

awk '{if (NR > 1) print '$tMIN'+(NR-1)*'$dt',$1}' $HFFILE |\
    gmt psxy $PROJ $LIMS -W1p -K -O >> $PSFILE
NPTAVG=100
NPTS=`wc $HFFILE | awk '{print $1-1}'`
awk '{
    if (NR>1) {
        hf[NR-1] = $1
    }
}END{
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
    gmt psxy $PROJ $LIMS -W1p,green -K -O >> $PSFILE



# Depth vs. temperature
echo "$0: plotting depth vs. temperature"
LIMS="-R$TMIN/$TMAX/$ZMIN/$ZMAX0"
PROJ="-JX5i/5i"

echo 0 0 | gmt psxy $PROJ $LIMS -X5.5i -Y-4i -K -O >> $PSFILE

gmt psbasemap $PROJ $LIMS -Bxa20g20+l"Temperature (\260C)" -Bya10g5+l"Depth (km)" -BWeSn -K -O >> $PSFILE

if [ "$GEOTHERM_FILE" != "" ]
then
    for T in $GEOTHERM_PLOT_TIME
    do
        T2=$(echo $T $tMIN | awk '{print $1+$2}')
        echo "Plotting geotherm at $T Ma since start of model, $T2 Ma until end of model"
        awk '{
            if (/>/) {
                p = 0
                if ($3==0 && $4=='$T') {
                    p = 1
                    print "> -W1p"
                } else if ($4=='$T') {
                    p = 1
                    print "> -W0.5p,105"
                }
            } else if (p==1) {
                print $1,-$2
            }
        }' $GEOTHERM_FILE > geotherm_$T.tmp
        WC_GEOTHERM=$(wc geotherm_$T.tmp | awk '{print $1}')
        if [ $WC_GEOTHERM -le 0 ]
        then
            echo "    Could not find a geotherm output at the specified time ($T Ma)!"
            echo "    Geotherms are available at the following times:"
            grep ">" $GEOTHERM_FILE | awk '{print "    " $4,"Ma"}'
        else
            gmt psxy geotherm_$T.tmp $PROJ $LIMS -K -O >> $PSFILE
            awk '{if(NR>1&&($2<='$ZMIN'||$1>='$TMAX')){print $1,$2,"8,2 LM '$T2' Ma";exit}}' geotherm_$T.tmp |\
                gmt pstext $PROJ $LIMS -F+f+j -D0.025i/0 -Gwhite -N -K -O >> $PSFILE
        fi
    done
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
echo 0 0 | gmt psxy $PROJ $LIMS -K -O -Y-4i >> $PSFILE
if [ "$GEOTHERM_FILE" != "" ]
then
    echo "$0: plotting temperature contours over time"
    LIMS="-R$tMIN/$tMAX/$ZMIN/$ZMAX"
    PROJ="-JX4i/3i -P"

    if [ "$TIMING_FILE" != "" ]; then plot_tectonic_timing; fi

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


# Closure temperature timing depth vs. time
echo 0 0 | gmt psxy $PROJ $LIMS -K -O -X5i >> $PSFILE
if [ "$CLOSURE_FILE" != "" ]
then
    echo "$0: plotting closure temperature timing"
    ZMIN2=$(awk '{if(!/>/){print $2}}' $CLOSURE_FILE | gmt gmtinfo -C | awk '{print $1}')
    ZMAX2=$(awk '{if(!/>/){print $2}}' $CLOSURE_FILE | gmt gmtinfo -C | awk '{print $2}')
    # tMIN2=$(awk '{if(!/>/){print '$tMIN'+$1}}' $CLOSURE_FILE | gmt gmtinfo -C | awk '{print $1}')
    # tMAX2=$(awk '{if(!/>/){print '$tMIN'+$1}}' $CLOSURE_FILE | gmt gmtinfo -C | awk '{print $2}')
    if [ "$CLOSURE_LIMS" == "" ]
    then
        LIMS="-R$tMIN/$tMAX/$ZMIN2/$ZMAX2"
    else
        LIMS=$CLOSURE_LIMS
    fi
    PROJ="-JX4i/3i -P"

    # if [ "$TIMING_FILE" != "" ]; then plot_tectonic_timing; fi

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
    gmt psbasemap $PROJ $LIMS -Bxa${AGE_TIKS}g${AGE_TIKS}+l"Age (Ma)" -Bya${DEP_TIKS}g${DEP_TIKS}+l"Depth (km)" -BWeSn -K -O >> $PSFILE

    awk 'BEGIN{p=0}{
        if ($1~/>/) {
            if (!/Partial/) {
                p = 1
                print $0
            } else {
                p = 0
            }
            getline
        }
        if (p==1) {
            if ($1!="#" && $1>'$dt') {
                print '$tMIN'+$1,$2
            }
        }
    }' $CLOSURE_FILE > closure.tmp
    gmt psxy closure.tmp $PROJ $LIMS -W1p -K -O >> $PSFILE
    gmt psxy closure.tmp $PROJ $LIMS -Sc2p -W1p -G155 -K -O >> $PSFILE

    awk 'BEGIN{p=0}{
        if ($1~/>/) {
            if (/Partial/) {
                p = 1
                print $0
            } else {
                p = 0
            }
            getline
        }
        if (p==1) {
            if ($1!="#" && $1>'$dt') {
                print '$tMIN'+$1,$2
            }
        }
    }' $CLOSURE_FILE > closure.tmp
    gmt psxy closure.tmp $PROJ $LIMS -G255/235/230 -K -O >> $PSFILE

    # Label curves
    awk 'BEGIN{p=0}{
        if (/>/) {
            temp = $2+0
            p = 1
            getline
        }
        if (!/#/ && $1>'$dt' && p==1) {
            print '$tMIN'+$1,$2,"10,2 CM",temp "\\260"
            p = 0
        }
    }' $CLOSURE_FILE |\
        gmt pstext $PROJ $LIMS -Gwhite -F+f+j -N -K -O >> $PSFILE

    awk '{if(/>/){print $0}else{print '$tMIN'+$1,$2}}' MWX_test6.3.trange | gmt psxy $PROJ $LIMS -W2p,blue -K -O >> $PSFILE

fi


echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

gmt psconvert $PSFILE -Tg -A
gmt psconvert $PSFILE -Tf -A
echo "Created file $(basename $PSFILE .ps).png"
echo "Created file $(basename $PSFILE .ps).pdf"

rm $PSFILE
rm gmt.*
rm color_list.tmp
