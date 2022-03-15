#!/bin/bash

# from https://gist.github.com/e7d/e43e6586c1c2ecb67ae2

export AnException=100
export AnotherException=101

source try-catch.sh

# start with a try
try
(   # open a subshell !!!
    echo "do something"
    [ someErrorCondition ] && throw $AnException

    echo "do something more"
    executeCommandThatMightFail || throw $AnotherException

    throwErrors # automaticatly end the try block, if command-result is non-null
    echo "now on to something completely different"
    executeCommandThatMightFail

    echo "it's a wonder we came so far"
    executeCommandThatFailsForSure || true # ignore a single failing command

    ignoreErrors # ignore failures of commands until further notice
    executeCommand1ThatFailsForSure
    local result = $(executeCommand2ThatFailsForSure)
    [ result != "expected error" ] && throw $AnException # ok, if it's not an expected error, we want to bail out!
    executeCommand3ThatFailsForSure

    echo "finished"
)
# directly after closing the subshell you need to connect a group to the catch using ||
catch || {
    # now you can handle
    case $ex_code in
        $AnException)
            echo "AnException was thrown"
        ;;
        $AnotherException)
            echo "AnotherException was thrown"
        ;;
        *)
            echo "An unexpected exception was thrown"
            throw $ex_code # you can rethrow the "exception" causing the script to exit if not caught
        ;;
    esac
}
