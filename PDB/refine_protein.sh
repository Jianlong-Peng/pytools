#!/bin/bash

if [ $# -ne 2 ];then
	echo "Usage: $0 start end"
	echo "  OBJ: to run 'refinestruct' in ./pdb{start} ... ./pdb{end}"
	exit
fi

for i in $(seq $1 $2)
do
	cd pdb${i}
	for inp_name in `ls *.inp`
	do
		echo "inp file: ./pdb${i}/${inp_name}"
		~/schrodinger/refinestruct ${inp_name} -WAIT
	done
	cd ..
done

