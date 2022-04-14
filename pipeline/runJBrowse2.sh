#!/bin/bash -l
#SBATCH -p nolimit -n 1

set -eu

# source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
	Usage: runJBrowse2.sh <dir> <port>
"

# usage when we singularity it
# "
# 	Usage: runJBrowse2.sh <singularity container> <dir> <port>
# "

# safety
# [[ $# -ne 3 ]] && abort "This script expects 3 arguments"
# [[ ! -f $1 ]] && abort "The first argument needs to be an existing singularity file"
# [[ ! -d $2 ]] && abort "The second argument needs to be an existing directory"
# [[ ! -f $2/config.json ]] && abort "The second argument needs to be an initialized jbrowse2 directory"
# [[ $3 -le 20000 ]] && abort "The port should be above or equal 20000"
# [[ $3 -gt 30000 ]] && abort "The port should be below 30000"

# safety
[[ $# -ne 2 ]] && abort "This script expects 2 arguments"
[[ ! -d $1 ]] && abort "The first argument needs to be an existing directory"
[[ ! -f $1/config.json ]] && abort "The first argument needs to be an initialized jbrowse2 directory"
[[ $2 -le 20000 ]] && abort "The port should be above or equal 20000"
[[ $2 -gt 30000 ]] && abort "The port should be below 30000"

# run
#docker run -d --rm -p $3:3000 -v $2:/var/www/html/jbrowse2 delhomme/upscb-jbrowse2
docker run -d --rm -p $2:3000 -v $1:/var/www/html/jbrowse2 delhomme/upscb-jbrowse2
