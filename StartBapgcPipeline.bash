#!/bin/bash

# Additional information:
# =======================
# Dit script is geschreven door Sebastiaan de Vriend en
# kan gebruikt worden bij de opdracht van bapgc door
# Jesse Kerkvliet.
# Het script is geschreven om de bapgc pipeline te starten.
# In het begin wordt er een temp en result directory gemaakt.
# Daarna wordt gevraagd welke pathway ge-analyseerd 
# kan worden. Het resultaat wordt geschreven in outputfiles en
# de tijdelijke bestanden worden opgeslagen in de directroy temp.
# Ook wordt er aan de user gevraagd of deze packages wilt
# Downloaden die gebruikt worden in de pipeline.
# Als de gebruiker geen packages heeft, wordt er aangeraden
# om hier ja op te beantwoorden.

if [ "$1" == "--h" ] || [ "$1" == "--help" ] || [ "$1" == "-h" ] || [ "$1" == "-help" ]
then
	echo "" 
	echo "Dit script is voor het opstarten van de bapgc pipeline."
	echo "Er wordt aan de gebruiker gevraagd welke pathway ge-analyseerd"
	echo "gaat worden en daarvoor wordt de BapgcPipeline mee gestart."
	echo "Daarnaast worden de directories temp en output gereset voor"
	echo "dataopslag."
	echo "Daarnaast is er de keuze op packages te downloaden."
	echo ""
	echo "Gebruikt als volgt:"
	echo ""
	echo "./StartBapgcPipeline.bash" 
	exit
fi

# Opties voor pathways.
# Pathways Jesse: hsa04722 hsa00980 hsa04920 hsa00360
opties="hsa04916 hsa00500 hsa00480" # Vul hier je pathways in.
echo "Welke pathway wil je analyseren?"
select optie in ${opties}
	do
	# Deze line wordt maximaal 3, voor de verschillende pathways.
	if [ "${REPLY}" = '1' ] || [ "${REPLY}" = '2' ] || [ "${REPLY}" = '3' ]
		then
			break
		else
			echo "Ongeldige keuze opgegeven. Voer een 1 in voor zoeken, 
			    voor een 2 in voor de analyse en 3 om te stoppen."
			echo "Voer je keuze in: "
			echo ${opties}
	fi
	done
boolPW=${optie}
boolJN="j n"
echo "Wil je de benodigde packages downloaden?"
select ant in ${boolJN}
    do
	if [ "${REPLY,,}" = 'j' ] || [ "${REPLY,,}" = 'n' ]
	
	then
	    boolPackages="${REPLY,,}"
		break
	fi	
    done
# Declaratie folder bestanden
outdir="$(pwd | awk '{ print $0 "/outputfiles"}')"${PW}
tempdir="$(pwd | awk '{ print $0 "/temp"}')"${PW}
#create new directories
mkdir ${outdir}
mkdir ${tempdir}
# Starten pipeline. Verander directory om lokale scripts aan te
# roepen.
cd Scripts/
sh ./BAPGC-Pipeline.sh ${outdir} ${tempdir} ${boolPackages} ${boolPW}
