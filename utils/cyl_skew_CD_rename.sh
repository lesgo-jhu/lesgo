#!/usr/bin/bash
FBASE="cylinder_skew_CD.dat"
FCASE="$1";#"ngen2.200k"
NPROC=$2;#128;
NGEN=$3;#2;
echo $FBASE
for (( np=0; np<=$NPROC-1; np++ ))
do
	for (( ng=1; ng <= $NGEN; ng++ ))
	do
		FTOT="$FBASE.g$ng.c$np";
		if [ -e "$FTOT" ]; then
			cp -v $FTOT $FTOT.$FCASE
		fi
	done
done


