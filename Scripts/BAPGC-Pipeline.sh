#!/bin/bash
if [ "$1" = "--h" ] || [ "$1" = "--help" ] || [ "$1" = "-h" ] || [ "$1" = "-help" ]
then
    echo ""
	echo "DOCUMENTATIE GEHAALD."
	echo "" 
	
    exit
fi
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
echo "Ophalen pathway informatie is klaar."
echo "Starten promoter script."
./PromoterSeqsFromGenes.R ${TempDir}
echo "Klaar met script"
