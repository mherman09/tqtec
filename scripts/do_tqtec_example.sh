#!/bin/bash

#####
#   PARSE INPUT
#####
BASENAME=$1
if [ "$BASENAME" == "" ]
then
    echo "$0: no model basename specified" 1>&2
    echo "Usage: $0 BASENAME" >&2
    exit 1
fi

if [ ! -f $BASENAME.in ]
then
    echo "$0: no input file found named $BASENAME.in" 1>&2
    echo "Usage: $0 BASENAME" >&2
    exit 1
fi


#####
#   RUN TQTEC
#####
echo "$0: Running tqtec"
../build/tqtec \
    -f $BASENAME.in \
    -o $BASENAME.out \
    -timing timing.tmp \
    -geotherm geotherm.tmp \
    -v 2 || exit 1


#####
#   RUN READTQTEC
#####
echo "$0: Running readtqtec"
../build/readtqtec \
    $BASENAME.out \
    -dep dep.tmp \
    -temp temp.tmp \
    -time time.tmp \
    -hf hf.tmp \
    -closure closure.tmp 60 120 || exit 1



#####
#   PLOT MODEL RESULTS
#####
# Horizon depth versus time
../scripts/plot_tqtec_dep_vs_time.sh dep.tmp \
    -timing timing.tmp
# Horizon temperature versus time
../scripts/plot_tqtec_temp_vs_time.sh temp.tmp \
    -timing timing.tmp
# Surface heat flow versus time
../scripts/plot_tqtec_hf_vs_time.sh hf.tmp \
    -timing timing.tmp
# Horizon temperature versus depth
../scripts/plot_tqtec_temp_vs_dep.sh dep.tmp temp.tmp \
    -geotherm geotherm.tmp \
    -geotherm:time 10 \
    -geotherm:time 20
# Horizon temperature contours superimposed on depth versus time
../scripts/plot_tqtec_temp_contours.sh dep.tmp geotherm.tmp \
    -timing timing.tmp
# Horizon depth versus time
../scripts/plot_tqtec_all.sh temp.tmp dep.tmp hf.tmp \
    -timing timing.tmp \
    -closure closure.tmp \
    -geotherm geotherm.tmp \
    -geotherm:time 10 \
    -geotherm:time 20 || exit

#####
#   CLEAN UP
#####
test -d ${BASENAME}_results || mkdir ${BASENAME}_results
cp *.pdf ${BASENAME}_results
mv $BASENAME.out ${BASENAME}_results
rm *.tmp
