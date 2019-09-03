module load bioinfo-tools
module load busco

# $0 [option] <name> <fasta> <out>
# option should be -m OGS or -m trans
# only runs on plants
MODE="trans"
cd $out
python3 $BUSCO_PATH/BUSCO_plants.py -o $1 -in $2 -l $BUSCO_DATA/plantae -m $MODE