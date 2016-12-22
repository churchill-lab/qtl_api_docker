#!/bin/bash

good=false
num_tries=0
RSCRIPT=''
RDATA=''

while [ "$good" = "false" ] && [ $num_tries -lt 50 ] ;
do
    echo "($num_tries) Trying...."

    if [ -d /api/data ] && [ -d /api/R ]; 
    then
        # the directory exists
        RDATA=`find /api/data -type f -name "*.R*ata"`
        RSCRIPT=`find /api/R -type f -name "*.R"`
        echo "Using R script: $RSCRIPT"
        echo "Using R data: $RDATA"

        if [ -n "$RDATA" ] && [ -n "$RSCRIPT" ];
        then
            echo "Files found"
            good=true
        fi
    else
        echo 'Paths do not exist'
    fi

    sleep 1
    num_tries=$((num_tries+1))
done

if [ "$good" = "true" ];
then
    echo "Using R script: $RSCRIPT"
    echo "Using R data: $RDATA"
    /usr/local/bin/Rscript /api/run_plumber.R $RSCRIPT $RDATA
fi

