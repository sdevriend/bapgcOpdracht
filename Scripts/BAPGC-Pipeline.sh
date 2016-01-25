#!/bin/bash

ResultDir=$1
TempDir=$2
PackageBool=$3
Pathway=$4

if [ "${PackageBool}" = 'j' ]
then
    ./DownloadRPackages.R
    echo "Klaar met downloaden
fi

echo "Klaar met script"
