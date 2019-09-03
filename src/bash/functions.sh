#!/bin/bash

# This file is meant to contain only functions to be sourced as
# source $UPSCb/src/bash/functions.sh

### ---------------------------------------------------------------------------
## logic functions

usage(){
  echo >&2 "$USAGETXT"
  exit 1;
}

abort(){
  echo >&2 $1
  usage
}

### ---------------------------------------------------------------------------
## array functions

# from https://stackoverflow.com/questions/3685970/check-if-a-bash-array-contains-a-value
containsElement () {
  local e
  for e in "${@:2}"; do
    [[ "$e" == "$1" ]] && echo 0 && return 0; 
  done
  echo 1
}

### ---------------------------------------------------------------------------
## preflight functions
isExec () {
  tool=`which $1 2>/dev/null`
  if [ ! -z $tool ] && [ -f $tool ] && [ -x $tool ]; then
    return 0
  else
    abort "The tool $tool is not available."
  fi
}

isEnvVarSet () {
  if [ ! -z $1 ]; then
    abort "The environment variable $1 is not set"
  fi
}

