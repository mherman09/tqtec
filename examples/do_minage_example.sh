#!/bin/bash


#####
#   USAGE
#####
function usage() {
    echo "Usage: $0 BASENAME -aft|-ahe"
    exit 1
}


#####
#   PARSE INPUT
#####
BASENAME=$1
if [ "$BASENAME" == "" ]
then
    echo "$0: no model basename specified" 1>&2
    usage
fi

if [ ! -f $BASENAME.in ]
then
    echo "$0: no input file found named $BASENAME.in" 1>&2
    usage
fi

MODE=$2
if [ "$MODE" == "-aft" ]
then
    echo "$0: modeling apatite fission tracks"
elif [ "$MODE" == "-ahe" ]
then
    echo "$0: modeling apatite (U-Th)/He"
else
    echo "$0: no minage mode named \"$MODE\"" 1>&2
    usage
fi





#####
# 
#####
if [ "$MODE" == "-aft" ]
then
    ../build/minage \
        -temp $BASENAME.in \
        -aft aft.tmp || exit 1
    ../scripts/plot_minage_aft_final.sh aft.tmp -temp $BASENAME.in
fi