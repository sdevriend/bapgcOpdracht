#!/bin/bash

# Additional information:
# =======================
# Dit script is geschreven door Sebastiaan de Vriend en
# kan gebruikt worden bij de opdracht van bapgc door
# Jesse Kerkvliet.
# Dit script roept losse script aan voor de BAPGC pipeline
# en zorgt ervoor dat alle tusenbestanden op correcte manier beheerd
# worden.  Het script kan aangeroepen worden met:
# ./BAPGC-Pipeline.sh ResultaatDirectory TijdelijkeDirectory j/n Pathwaynaam.
# 
# Het script loopt alle tussenstappen af voor het maken van een verslag
# en de daarbij behorende data te genereren met diverse
# scripts. Op het laatst wordt er een phylogenetische boom
# gemaakt met ClustalW door fasta bestanden te mergen in een enkele
# fasta file.

if [ "$1" = "--h" ] || [ "$1" = "--help" ] || [ "$1" = "-h" ] || [ "$1" = "-help" ]
then
    echo ""
	echo " Dit script roept losse script aan voor de BAPGC pipeline"
    echo " en zorgt ervoor dat alle tusenbestanden op correcte manier "
	echo "beheerd worden. Het script kan aangeroepen worden met:"
	echo "./BAPGC-Pipeline.sh ResultaatDirectory TijdelijkeDirectory j/n Pathwaynaam."
	echo ""
	echo "Het script maakt op het einde een phylogenetische boom"
	echo " doormiddel van clustalW."
    exit
fi
# Declaratie bestanden en namen.
ResultDir=$1
TempDir=$2
PackageBool=$3
Pathway=$4
# Optie om alle packages te downloaden. Hiervoor zijn root rechten nodig.
if [ "${PackageBool}" = 'j' ]
then
    echo "Root permissie is nodig voor package installatie."
    sudo ./DownloadRPackages.R
    echo "Klaar met downloaden"
fi
echo "Ophalen Pathway info van Kegg."
#./GenesFromPathway.r ${TempDir} ${Pathway}
echo "Ophalen pathway informatie is klaar."
echo "Starten promoter script."
#./PromoterSeqsFromGenes.R ${TempDir}
echo "Promotor sequenties zijn gevonden. Door met motifs zoeken."
#./MotifsFromPromoterSequences.R ${TempDir}
echo "Genen zoeken..."
#./GenesFromMotifs.R ${TempDir}
echo "Genen gevonden. Zoeken naar genen met functie."
#./EnzymesWithFunctions.R ${TempDir}
echo "Enzymen gevonden. Bepalen Hydrofobic eigenschappen..."
#./Hydrofobiciteit.R ${TempDir}
echo "Hydrofobic eigenschappen opgezocht."
echo "Door met intron/exon bepaling."
#./IntronsAndExons.R ${TempDir}
echo "Maken MSA."
#./MSA.R ${TempDir}
echo "Testen phylo"
#./PhylogenyOfParalogs.R ${TempDir}
# Maken van multi fasta met de volgende twee pipes:
#cat ${TempDir}/*.fasta | awk '{if (substr($1,1,1) == ">") print $0" "; else print $1}' | tr -d "\n" | tr ">" "\n"  > ${TempDir}/multitemp.fa
#cat ${TempDir}/multitemp.fa  | awk '{if(NR>1)print $0}' | awk '{print ">"$0}' | sed s/" "/"\n"/  > ${TempDir}/tree.fa
#clustalw data -infile=${TempDir}/tree.fa -tree -outputtree=dist 
#./generateTree.R ${TempDir}
echo "Klaar met script. Door met verslag."
cp ${TempDir}/Hydrofobiciteit.png .
cp ${TempDir}/IntronsExons.png .
cp ${TempDir}/MSARESULT.txt .
cp ${TempDir}/Pathway.png .
cp ${TempDir}/tree.png .
cp ${TempDir}/VennDiagram.png .
cp ${TempDir}/FoundMotifs.csv .
# Run latex script

mv Hydrofobiciteit.png ${ResultDir}
mv IntronsExons.png ${ResultDir}
mv MSARESULT.txt ${ResultDir}
mv Pathway.png ${ResultDir}
mv tree.png ${ResultDir}
mv VennDiagram.png ${ResultDir}
mv FoundMotifs.csv ${ResultDir}
