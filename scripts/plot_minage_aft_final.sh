#!/bin/bash

#####
#   USAGE STATEMENT
#####
function usage() {
    echo "Usage: $0 AFTFILE [-temp TEMPFILE] [-timing TIMINGFILE]" 1>&2
    exit 1
}


#####
#   PARSE COMMAND LINE
#####
# First argument must be AFT file from minage
AFT_FILE=$1

# Check that file is defined and exists
if [ "$AFT_FILE" == "" ]
then
    echo "$0: no AFT output file defined" 1>&2
    usage
fi
if [ ! -f $AFT_FILE ]
then
    echo "$0: no AFT file found named \"$AFT_FILE\"" 1>&2
    usage
fi
shift

# Parse other command line arguments
TEMP_FILE=
TIMING_FILE=
while [ "$1" != "" ]
do
    case $1 in
        -temp) shift; TEMP_FILE=$1;;
        -timing) shift; TIMING_FILE=$1;;
        *) usage;;
    esac
    shift
done



#####
#   EXTRACT CONTENTS OF AFT FILE
#####

# Retention age
AFT_RETENTION_AGES=`awk '{
    if (/RETENTION AGE/) {
        getline
        print $0
        exit
    }
}' $AFT_FILE`

# Fission track age
AFT_AGES=`awk '{
    if (/FISSION TRACK AGE/) {
        getline
        print $0
        exit
    }
}' $AFT_FILE`

# Corrected length distribution
awk 'BEGIN{p=0}{
    if ($1 == "#") {
        if (/Final \(corrected\) track lengths/) {
            p = 1
            getline
        } else {
            p = 0
        }
    } 
    if (p==1) {
        print $0
    }
}' $AFT_FILE > aft_corrected_hist.tmp

# Max track count in a bin
NMAX=`awk '{for(i=2;i<=NF;i++){print $i}}' aft_corrected_hist.tmp | gmt gmtinfo -C | awk '{print $2}'`



#####
#   PLOT RESULTS
#####
gmt set PS_MEDIA 100ix100i

PSFILE=minage_aft_final.ps

PROJ=-JX6i/3i
CTMAX=`echo $NMAX | awk '{print $1*1.12}'`
LIMS=-R0/20/0/$CTMAX

# Colors for each horizon
NCOL=$(sed -ne "2p" aft_corrected_hist.tmp | awk '{print NF-1}')
NCOL2=$(echo $NCOL 1.1 | awk '{print $1*$2}')
gmt makecpt -Cmagma -T0/$NCOL2/1 > color_list.cpt
awk '{if(NF==5){print $2}}' color_list.cpt > color_list.tmp


# Initialize length distribution plot
gmt psxy -T -K -Y95i > $PSFILE


# Different panel for each distribution
for COL in `seq 1 $NCOL`
do

    # Set color for this horizon
    COLOR=$(sed -ne "${COL}p" color_list.tmp)

    # Plot gridlines
    gmt psbasemap $PROJ $LIMS -Byg -BSW -K -O \
        --MAP_GRID_PEN_PRIMARY=0.5p,205,4_2:0 >> $PSFILE

    # Mean and median lengths
    MEAN=`awk 'BEGIN{i='$COL'+1;sum=0;n=0}{
        n = n + $i
        len = $1
        sum = sum + $i*len
    }END{print sum/n}' aft_corrected_hist.tmp`
    echo $MEAN $CTMAX | awk '{print $1,0;print $1,$2}' |\
        gmt psxy $PROJ $LIMS -W2p,$COLOR -K -O >> $PSFILE
    MEDIAN=`awk 'BEGIN{i='$COL'+1;n=0}{
        for (j=1;j<=$i;j++) {
            n++
            len[n] = $1
        }
    }END{print len[int(n/2)]}' aft_corrected_hist.tmp`
    echo $MEDIAN $CTMAX | awk '{print $1,0;print $1,$2}' |\
        gmt psxy $PROJ $LIMS -W2p,$COLOR,2_6:0 -K -O --PS_LINE_CAP=round >> $PSFILE

    # Plot length distributions
    awk 'BEGIN{i='$COL'+1}{
        ct = $i
        print ">"
        print $1-0.5,0
        print $1-0.5,$i
        print $1+0.5,$i
        print $1+0.5,0
        print $1-0.5,0
    }' aft_corrected_hist.tmp |\
        gmt psxy $PROJ $LIMS -W1p -G$COLOR -K -O >> $PSFILE


    # Label length distributions with ages
    # Fission track age
    AFT_AGE=`echo $AFT_AGES | awk '{printf("%.1f"),$'$COL'}'`
    echo "Fission Track Age = $AFT_AGE Ma" |\
        gmt pstext $PROJ $LIMS -F+f12,3+cLB -D0.05i/19p -K -O >> $PSFILE

    # Retention age
    AFT_RETENTION_AGE=`echo $AFT_RETENTION_AGES | awk '{printf("%.1f"),$'$COL'}'`
    echo "Retention Age = $AFT_RETENTION_AGE Ma" |\
        gmt pstext $PROJ $LIMS -F+f12,3+cLB -D0.05i/5p -K -O >> $PSFILE

    # Plot & label frame
    gmt psbasemap $PROJ $LIMS -Bxa1+l"Track Length (microns)" -Byf+l"Count" -BWeSn -K -O >> $PSFILE
    echo Horizon $COL | gmt pstext $PROJ $LIMS -F+f20,3+cRT -D-0.05i/-0.04i -K -O >> $PSFILE

    #####
    #   PLOT TEMPERATURE HISTORY
    #####
    if [ "$TEMP_FILE" != "" ]
    then
        awk '{
            if(NR==1){
                dt=$2
            }else{
                for(i=1;i<='$NCOL';i++){
                    print (NR-1)*dt,$i
                }
            }
        }' $TEMP_FILE | gmt gmtinfo -C > minmax.tmp
        TIMEMIN=`awk '{print $1}' minmax.tmp`
        TIMEMAX=`awk '{print $2}' minmax.tmp`
        TEMPMIN=`awk '{print $3-5}' minmax.tmp`
        TEMPMAX=`awk '{print $4+5}' minmax.tmp`
        PROJ_TEMP=-JX2.0i/1.6i
        PROJ_TEMP_AXIS=-JX-2.0i/1.6i
        LIMS_TEMP=-R0/$TIMEMAX/$TEMPMIN/$TEMPMAX
        SHFT="-Xa0.2i -Ya1.2i"

        # White background
        echo $LIMS_TEMP | sed -e "s/-R//" |\
            awk -F/ '{print $1,$3;print $1,$4;print $2,$4;print $2,$3}' |\
            gmt psxy $PROJ_TEMP $LIMS_TEMP -Gwhite -K -O $SHFT >> $PSFILE

        # Tectonic action timing
        YMIN=$(echo $LIMS_TEMP | sed -e "s/-R//" | awk -F/ '{print $3}')
        YMAX=$(echo $LIMS_TEMP | sed -e "s/-R//" | awk -F/ '{print $4}')
        awk 'BEGIN{print ">"}{
            if ($1=="burial" || $1=="uplift" || $1=="thicken") {
                print ">"
                print 0+$3,'$YMIN'
                print 0+$4,'$YMIN'
                print 0+$4,'$YMAX'
                print 0+$3,'$YMAX'
                print 0+$3,'$YMIN'
            }
        }' $TIMING_FILE |\
            gmt psxy $PROJ_TEMP $LIMS_TEMP -G245 -K -O $SHFT >> $PSFILE
        awk 'BEGIN{print ">"}{
            if ($1=="thrust") {
                print ">"
                print 0+$3,'$YMIN'
                print 0+$3,'$YMAX'
            }
        }' $TIMING_FILE |\
            gmt psxy $PROJ_TEMP $LIMS_TEMP -W3p,225 -K -O $SHFT >> $PSFILE

        # Temperature versus time
        awk '{
            if (NR==1) {
                dt = $2
            } else {
                print (NR-1)*dt,$'$COL'
            }
        }' $TEMP_FILE |\
            gmt psxy $PROJ_TEMP $LIMS_TEMP -W1p,$COLOR -K -O $SHFT >> $PSFILE
        gmt psbasemap $PROJ_TEMP_AXIS $LIMS_TEMP -BwESn \
            -Bxa+l"Time Before Present (Ma)" -Bya+l"Temperature (\260C)" -K -O $SHFT \
            --FONT_ANNOT_PRIMARY=8 --FONT_LABEL=10 \
            --MAP_LABEL_OFFSET=5p --MAP_ANNOT_OFFSET_PRIMARY=3p --MAP_TICK_LENGTH=3p >> $PSFILE
    fi


    gmt psxy -T -K -O -Y-4.0i >> $PSFILE
done




gmt psxy -T -O >> $PSFILE
gmt psconvert $PSFILE -Tf -A


#####
# CLEAN UP
#####
rm $PSFILE
rm *.tmp
rm color_list.cpt