#!/bin/bash
#SBATCH -p core
#SBATCH -n 1

set -ex

usage(){
  echo >&2 \
  "
  Usage: $(basename $0) <function> [options] -- docker args
  
  Functions: 
    run:   run a docker image, requires the -i option and any additional arguments
    list:  list all running containters, better used as a 'srun' command
    clean: stop and remove a docker container, requires the -k option (the docker 
           container id, which can be retrieved using the 'list' function above)
  
  Options:
            -c a command to execute if not the image default
            -i the image name
            -k the container id
            -d dry-run simply print the command
            
  Note: The '--' argument is a special argument that forward remaining arguments to the docker run command
        DO NOT PROVIDE THE docker run -d option as an extra command line, that would put the SLURM queue out
        of sync with the available resources
  "
  exit 1
}

# exit if user not in docker group
ALWD=`groups | grep -c docker`

if [ $ALWD -ne 1 ]; then
  echo "Your user is not part of the docker group, contact your administrator."
  usage
fi

## VARS
DRY=0
IMAGE=
CID=
CMD=

## get the options
while getopts c:di:k: option
do
  case "$option" in
      c) CMD=$OPTARG;;
      d) DRY=1;;
	    i) IMAGE=$OPTARG;;
	    k) CID=$OPTARG;;
	    \?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## arguments
case "$1" in
    run)
    
    # get the cmd
    #CMD=$1
    shift
    
    # get additional options
    if [ $# != 0 ]; then
  	## drop the --
	    shift
    fi

    # check the image
    if [ -z $IMAGE ]; then
      echo "You need to provide the -i argument to the run command"
    fi

    EXEC="docker run $@ $IMAGE $CMD"
    ;;
    list)
    
    # get the cmd
    CMD=$1
    shift
    
    # check for additional options
    if [ $# != 0 ]; then
      echo "You cannot give additional arguments to the 'list' command"
	    usage
    fi
    
    EXEC="docker ps -a"
    
    ;;
    clean)
    
    # get the cmd
    #CMD=$1
    shift
    
    # check for additional options
    if [ $# != 0 ]; then
      echo "You cannot give additional arguments to the 'list' command"
	    usage
    fi
    
    # check the container id
    if [ -z $CID ]; then
      echo "You need to provide the -k argument to the clean command"
    fi

    EXEC="docker stop $CID && docker rm $CID"
    
    ;;
    *)
	echo "Unsupported command. RTFM."
	usage;;
esac

if [ $DRY -eq 1 ]; then
  echo $EXEC
else
  exec $EXEC
fi
