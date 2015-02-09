#! /bin/bash

exit_code=0;

if [ $# -ne 0 ] && [ "$1" == "clean" ]
then
    rm -rf run/*.exe bin/*.o bin/*.a bin/*.d *.exe *.out
    ./run/remove_backups.sh
    exit_code=$?
else
    tmp_file=mktemp
    
    make -j 4 -k -r -R 2> >(tee $tmp_file >&2)
    exit_code=$?
    
    echo
    
    if [[ $exit_code != 0 ]] ; then
	echo "ERRORS AND WARNINGS:"
	cat $tmp_file >&2
    else
	echo "Compiled successfully without errors or warnings!"
    fi
    
    rm -rf $tmp_file
fi

exit $exit_code
