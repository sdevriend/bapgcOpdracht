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


#create new directories

outdir="$(pwd | awk '{ print $0 "/outputfiles"}')"
tempdir="$(pwd | awk '{ print $0 "/temp"}')"
opties="hsa04916 hsa00500 hsa00480 hsa04722 hsa00980 hsa04920 hsa00360" # Vul hier je pathways in.
echo "Welke pathway wil je analyseren?"
select optie in ${opties}
	do
		if [ "${REPLY}" = '1' ] || [ "${REPLY}" = '2' ] || [ "${REPLY}" = '3' ] || [ "${REPLY}" = '4' ] || [ "${REPLY}" = '5' ] || [ "${REPLY}" = '6' ]
		then
			break
		else
			echo "Ongeldige keuze opgegeven. Voer een 1 in voor zoeken, voor een 2 in voor de analyse en 3 om te stoppen."
			echo "Voer je keuze in: "
			echo ${opties}
		fi
	done
PW=${optie}
JN="j n"
echo "Wil je de benodigde packages downloaden?"
select ant in ${JN}
    do
	if [ "${REPLY,,}" = 'j' ] || [ "${REPLY,,}" = 'n' ]
	
	then
	    Packages="${REPLY,,}"
		break
	fi
		
    done
rm -rf outputfiles${PW}
rm -rf temp${PW}
mkdir outputfiles${PW}

mkdir temp${PW}

sh ./Scripts/BAPGC-Pipeline.sh ${outdir} ${tempdir} ${Packages} ${PW}
