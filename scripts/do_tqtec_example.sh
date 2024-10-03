#!/bin/bash

#####
#   USAGE STATEMENT
#####
function usage() {
    echo "Usage: $0 MODEL_FILE [...options...]" 1>&2
    echo "MODEL_FILE         Input TQTec model file" 1>&2
    echo "-plot:all          Generate a single figure file with all graphs" 1>&2
    echo "-plot:individual   Generate an individual figure file for each graph" 1>&2
    exit
}


#####
#   PARSE INPUTS
#####
# Check for input model file argument
MODEL_FILE=$1
if [ "$MODEL_FILE" == "" ]
then
    echo "$0: no model input file specified" 1>&2
    usage
fi

# Check whether input model file exists
if [ ! -f $MODEL_FILE ]
then
    echo "$0: no input file found named $MODEL_FILE" 1>&2
    usage
fi

# Parse arguments
shift
PLOT=ALL
while [ "$1" != "" ]
do
    case "$1" in
        -plot:all) PLOT="ALL";;
        -plot:individual) PLOT="INDIVIDUAL";;
        *) echo "No option $1" 1>&2; usage;;
    esac
    shift
done


# Extract model basename
MODEL_BASENAME=`basename $MODEL_FILE .in`


#####
#   RUN TQTEC
#####
echo "$0: Running tqtec"
tqtec \
    -f $MODEL_FILE \
    -o $MODEL_BASENAME.out \
    -timing timing.tmp \
    -geotherm geotherm.tmp \
    -v 2 || exit 1


#####
#   RUN READTQTEC
#####
echo "$0: Running readtqtec"
readtqtec \
    $MODEL_BASENAME.out \
    -dep dep.tmp \
    -temp temp.tmp \
    -time time.tmp \
    -hf hf.tmp \
    -closure closure.tmp 60 120 || exit 1



#####
#   PLOT MODEL RESULTS
#####
# Determine times to plot geotherms
MODEL_TIME=`head -1 temp.tmp | awk '{print -$1}'`
GEOTHERM_TIME_FLAG=
for T in `seq 10 10 $MODEL_TIME`
do
    GEOTHERM_TIME_FLAG="$GEOTHERM_TIME_FLAG -geotherm:time $T"
done

# Plot individual figure files or single figure file
if [ "$PLOT" == "INDIVIUDAL" ]
then
    # Horizon depth versus time
    plot_tqtec_dep_vs_time.sh dep.tmp \
        -timing timing.tmp
    # Horizon temperature versus time
    plot_tqtec_temp_vs_time.sh temp.tmp \
        -timing timing.tmp
    # Surface heat flow versus time
    plot_tqtec_hf_vs_time.sh hf.tmp \
        -timing timing.tmp
    # Horizon temperature versus depth
    plot_tqtec_temp_vs_dep.sh dep.tmp temp.tmp \
        -geotherm geotherm.tmp \
        $GEOTHERM_TIME_FLAG
    # Horizon temperature contours superimposed on depth versus time
    plot_tqtec_temp_contours.sh dep.tmp geotherm.tmp \
        -timing timing.tmp
elif [ "$PLOT" == "ALL" ]
then
    # Horizon depth versus time
    plot_tqtec_all.sh temp.tmp dep.tmp hf.tmp \
        -timing timing.tmp \
        -closure closure.tmp \
        -geotherm geotherm.tmp \
        $GEOTHERM_TIME_FLAG || exit
fi


#####
#   CLEAN UP
#####
test -d ${MODEL_BASENAME}_results || mkdir ${MODEL_BASENAME}_results
cp *.pdf ${MODEL_BASENAME}_results
mv $MODEL_BASENAME.out ${MODEL_BASENAME}_results
rm *.tmp
