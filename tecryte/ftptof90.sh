#!/bin/bash
FPREFIX=("write_real_data" "write_real_data_1d" "write_real_data_2d" "write_real_data_3d");
DATATYPE=("single" "double");
DATAPREC=("4" "8");
DATAFORMAT=( "data_format_single" "data_format_double" );

#if [ $# -ne 1 ]
#then
#  echo "Usage: `basename $0` {prefix}"
#  exit 1
#fi

# OS check
OSNAME=`uname`

# OS specific work arounds
if [[ "$OSNAME" == 'Darwin' || "$OSNAME" == 'FreeBSD' ]]; then
    SED_CMD='sed -i .tmp'
else
    SED_CMD='sed -i.tmp'
fi


for (( j = 0 ; j < ${#FPREFIX[@]} ; j++ ))
do 
  F90FILE=${FPREFIX[$j]}.f90
  FTPFILE=${FPREFIX[$j]}.ftp

  if [ ! -f $FTPFILE ]
  then
      echo template file $FTPFILE does not exist
      exit 2
  fi

  rm -f $F90FILE

  for (( i = 0 ; i < ${#DATATYPE[@]} ; i++ ))
  do
    cat  $FTPFILE >> $F90FILE
    $SED_CMD 's/<DATATYPE>/'${DATATYPE[$i]}'/g' $F90FILE
    $SED_CMD 's/<DATAPREC>/'${DATAPREC[$i]}'/g' $F90FILE
    $SED_CMD 's/<DATAFORMAT>/'${DATAFORMAT[$i]}'/g' $F90FILE
  done
  # Clean up
  rm -f $F90FILE.tmp
done



