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
#   PLOT TECTONIC-THERMAL MODEL RESULTS
#####
# All plots
../scripts/plot_tqtec_all.sh temp.tmp dep.tmp hf.tmp \
    -timing timing.tmp \
    -closure closure.tmp \
    -geotherm geotherm.tmp \
    -geotherm:time 10 \
    -geotherm:time 20 \
    -geotherm:time 30 \
    -geotherm:time 40 \
    -geotherm:time 50 || exit



#####
#   PLOT APATITE FISSION TRACK RESULTS
#####
../build/minage \
    -temp temp.tmp \
    -dep dep.tmp \
    -aft aft.tmp || exit 1
../scripts/plot_minage_aft_final.sh aft.tmp \
    -temp temp.tmp \
    -timing timing.tmp


#####
#   CLEAN UP
#####
test -d ${BASENAME}_results || mkdir ${BASENAME}_results
cp *.pdf ${BASENAME}_results
mv $BASENAME.out ${BASENAME}_results
rm *.tmp
