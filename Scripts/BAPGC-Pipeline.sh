#!/bin/bash

ResultDir=$1
TempDir=$2
PackageBool=$3
Pathway=$4
cd Scripts/
pwd
if [ "${PackageBool}" = 'j' ]
then
    echo "Root permissie is nodig voor package installatie."
    sudo ./DownloadRPackages.R
    echo "Klaar met downloaden"
fi
echo "Ophalen Pathway info van kegg."
./GenesFromPathway.r ${TempDir} ${Pathway}
echo "Klaar met script"
