#!/bin/bash

# Additional information:
# =======================
# Dit script is geschreven door Sebastiaan de Vriend en
# kan gebruikt worden bij de opdracht van bapgc door
# Jesse Kerkvliet.
# Het script is geschreven om de bapgc pipeline te starten.
# In het begin worden alle directories geleegd en daarna opnieuw
# aangemaakt. Daarna wordt gevraagd welke pathway ge-analyseerd 
# kan worden. Het resultaat wordt geschreven in outputfiles en
# de tijdelijke bestanden worden opgeslagen in de directroy temp.
# 

if [ "$1" == "--h" ] || [ "$1" == "--help" ] || [ "$1" == "-h" ] || [ "$1" == "-help" ]
then
	echo "" 
	echo "Dit script is voor het opstarten van de bapgc pipeline."
	echo "Er wordt aan de gebruiker gevraagd welke pathway ge-analyseerd"
	echo "gaat worden en daarvoor wordt de BapgcPipeline mee gestart."
	echo "Daarnaast worden de directories temp en output gereset voor"
	echo "dataopslag. "
	echo ""
	echo "Gebruikt als volgt:"
	echo ""
	echo "./StartBapgcPipeline.bash" 
	exit
fi

#clear directories
rm -rf ﻿outputfiles
rm -rf temp
sleep 1
#create new directories
mkdir outputfiles
sleep 1
mkdir temp
outdir="$(pwd | awk '{ print $0 "/outputfiles"}')"
outdir=
echo ${outdir}
opties="hsa04916 hsa00500"
select optie in ${opties}
	do
		if [ "${REPLY}" = '1' ] || [ "${REPLY}" = '2' ]
		then
			break
		else
			echo "Ongeldige keuze opgegeven. Voer een 1 in voor zoeken, voor een 2 in voor de analyse en 3 om te stoppen."
			echo "Voer je keuze in: "
			echo ${opties}
		fi
	done
sh ./BGSE-NGS-Pipeline.sh
