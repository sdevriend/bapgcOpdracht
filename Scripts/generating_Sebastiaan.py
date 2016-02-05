#!/usr/bin/python
# -*- coding: cp1252 -*-
import os
import sys
import argparse
import subprocess
import datetime

def showUsageInformation():
    print() 
    print("Het script bevat de strings nodig om een LaTeX bestand")
    print("Te genereren. Dit bestand wordt weggeschreven")
    print("")
    print("generating_Sebastiaan.py -n 'S1082434' -name 'Piet' -title 'Hello World'")
    print()
    sys.exit()

 

def front_page(argdict):
    '''
    De voorpagina wordt gegenereerd met het plaatje Pathway,
    de meegegeven argumenten (naam, studentnummer, titel) en
    de huidige maand.
    '''
    now = datetime.datetime.now()
    firstpage=r'''
    \documentclass{article}
    \usepackage{graphicx}
    \begin{document}
    \textbf{\Huge %(title)s \\}
    \vspace{1cm}
    \begin{figure}[h]
    \centering
    \includegraphics[height=10cm, width=10cm]{Pathway}
    \caption{Pathway}
    \end{figure}
    
    \vspace{3cm}
    \textbf{\large %(name)s \\}
    \vspace{1cm}
    \textbf{\normalsize %(student_num)s \\}
    \textbf{\normalsize %(month)s \\}
    '''
    month_year = "%s %s"%(str(now.strftime("%B")),str(now.year))
    argdict["month"] = month_year
    frontpage_Latex = firstpage%argdict
    return frontpage_Latex

def inleiding():
    '''
    De inleiding wordt gegenereerd. Dit is platte tekst en gebruikt geen
    variabelen
    '''
    inleiding=r'''
    \newpage
    {\Huge Inleiding}
    \newline
    \newline
    {Dit document wordt automatisch gegenereerd door de
    BapgcCoregulationPipeline. Het document gaat over de regulatie van
    genen in een metabolische route. Er wordt onderzocht welke genen
    een rol spelen in het Melanogenesis signaal-pathway.
    Melanogenesis is de productie van melanine. 
    Melanine is een pigment die gevonden wordt in de huid en het haar. 
    De taak van melanine is om de onderste laag van de huid te beschermen
    tegen schadelijke UV-B straling. 
    Melanine is donker opdat het meer licht absorbeert.
    Melanogenesis zorgt er dus voor dat de huid van een mens donkerder kleurt
    na blootstelling aan zonlicht. De onderzoeksvraag die gesteld wordt bij
    dit onderzoek is als volgt: Zijn er genen die hetzelfde gereguleerd
    worden als de genen van het Melanogenesis signaal-pathway?
    '''
    
    return inleiding

def materiaal_methode():
    '''
    De materiaal en methode wordt gegenereerd.
    Dit is platte tekst en gebruikt geen variabelen.
    '''
    materiaal_methode = r'''
        \newpage
        {\Huge Materiaal en Methode}
        \newline
        {Om de transcriptiefactoren van deze pathway te vinden, is voor elk van
        deze genen gezocht naar motifs die bij transcriptiefactoren horen.
        Hiervoor is een downstream regio van 3000 en en upstream regio van 300
        basenparen gebruikt. De gevonden motifs vormen samen de set
        transcriptiefactoren die de pathway reguleert.
        Hierna wordt het hele humane genoom gescand op deze motifs.
        De genen die elk van deze motifs bevatten, worden meegenomen als
        medegereguleerde genen. Dit zijn 19476 genen.
        }
        \newline
        \newline
        {Vijf biologische vragen zijn beantwoord met deze set genen.}
        \newline
        {1: Er wordt onderzocht hoe veel van de beste 100 medegereguleerde
        genen enzymen zijn en hoe veel van deze enzymen een bekende functie
        heeft.}
	\newline
        {2:	Er wordt een Multiple Sequence Alignment(MSA) uitgevoerd van
        de beste 50 medegereguleerde genen. Er wordt berekend welk percentage
        van het langste gen geconserveerd is in deze MSA.}
        \newline
        {3: 	Er wordt gekeken naar de hydrofobiciteit van de eiwitten van
        de beste 10 medegereguleerde genen. Dit is opgedeeld in de categorieën
        hydrofoob, hydrofiel en neutraal.}
        \newline
        {4: 	De vier meest verwante paralogen van de beste 10
        medegereguleerde genen worden gezocht. Hier wordt een fylogenetische
        analyse op toegepast.}
        \newline
        {5:	Er wordt gekeken naar de intron- en exonlengten van de beste 20
        medegereguleerde genen.}
        \newline
    '''
    return(materiaal_methode)

def resultaat(rMotifs):
    '''
    De functie genereert de resultatensectie. Als variabelen
    wordt een lijst met motifs en een gennaam + bijbehorend percentage
    meegegeven. De plaatjes worden uit de working directory gehaald
    door LateX
    '''
    result = r'''
    \newpage
    {\Huge Resultaten}
    {Uit het zoeken naar transcriptiefactoren zijn de volgende motifs
    geselecteerd:
    %s , %s, %s.}
    \newline
    {Dit resulteert in een set van 19476 genen}
    {Van de bovenste 100 genen van deze set is gezocht naar enzymen. 
    De hoeveelheid enzymen en de enzymen zonder functie is in figuur 2.}

    \begin{figure}[h]
    \includegraphics[width=0.5\textwidth]{VennDiagram}
    \caption{VennDiagram van alle enzymen in de dataset}
    \end{figure}

    {Van de bovenste 50 medegereguleerde genen is een MSA gemaakt.
    Hieruit is afgeleid dat
    gen %s voor %s geconserveerd is in deze MSA.
    Van de bovenste 10 medegereguleerde genen is het eiwit bekeken en berekend
    hoe veel aminozuren hydrofiel, hydrofoob danwel neutraal zijn.
    Deze percentages zijn te zien in
    figuur 3.}

    \begin{figure}[h]
    \includegraphics[width=0.5\textwidth]{Hydrofobiciteit}
    \caption{Hydrofobiciteit van de 20 beste genen.}
    \end{figure}

    {Van de bovenste 10 medegereguleerde genen is gezocht naar maximaal vier
    paralogen. Deze vier per gen zijn weergegeven in een fylogenetische boom,
    zoals te zien in figuur 4}

    \begin{figure}[h]
    \includegraphics[width=0.8\textwidth]{tree}
    \caption{{Phylogenetische boom van de 19 beste genen.}}
    \end{figure}

    {Van de bovenste 19 medegereguleerde genen is de lengte van de intronen en exonen bepaald.
    De resultaten hiervan zijn te zien in figuur 5}

    \begin{figure}[h]
    \includegraphics[width=0.8\textwidth]{IntronsExons}
    \caption{Intronen en exonen lengtes.}
    \end{figure}
    '''
    result = result % tuple(rMotifs)
    return(result)

def discussie_conclusie():
    '''
    De discussie en conclusie (disclusie) sectie wordt gegenereerd.
    Hier wordt platte tekst voor gebruikt en er worden geen variabelen
    gebruikt.
    '''
    disclusie = r'''
    \newpage
    \newpage
    {\Huge Discussie en conclusie}
    \newline
    Wat meteen duidelijk te zien is in de resultaten, is dat er 19476 genen
    gereguleerd worden. Deze hoeveelheid is erg groot voor een specifieke
    pathway. Mogelijk komt dit omdat er niet goed gefilterd is bij het zoeken
    naar genen die bij motieven horen en is er niet gekeken naar de
    specifieke frameworks die horen bij het pathway. Daarnaast zijn er drie
    motifs gevonden voor alle genen op het pathway. Dit had aangepast kunnen
    worden naar tien motifs die 90% van de genen reguleren. Dit had waarschijnlijk
    miner genen opgeleverd. Hierdoor is er te concluderen dat de gevonden genen
    niet specifiek genoeg zijn voor de onderzoeksvraag.
    In figuur 2 is te zien dat er drie enzymen zijn gevonden in de genenlijst.
    De drie enzymen hebben een bekende functie. In figuur 3 is goed te zien
    dat de gevonden eiwitten vooral hydrofoob zijn. In figuur 4 is te zien dat
    er geen concluderend resultaat is voor de fylogenetische boom, maar wel is
    er te zien dat er geen out-groups zijn.
    In figuur 5 is er te zien dat er veel genen laag op de intron as staan
    en laag op de exon as. Dit betekend dat veel genen kort zijn. Er wordt
    namelijk verwacht dat de meeste genen hoog op de intron as staan en
    op gemiddelde hoogte op de exon as.
    De onderzoeksvraag: Zijn er genen die hetzelfde gereguleerd
    worden als de genen van het Melanogenesis signaal-pathway? kan niet goed 
    beantwoord worden.
    \end{document}
    '''
    return disclusie

def readCSV():
    '''
    Het bestand met de gevonden motifs wordt geopend
    en geformatteerd. Een lijst met motifs wordt teruggegeven.
    ''' 
    openFile = open('FoundMotifs.csv','r')
    lines = openFile.readlines()
    openFile.close()
    motiflist = []
    for x in range(len(lines)):
        motif = lines[x].split(",")[1].strip("\n")
        motif = motif.strip("\r")
        if x > 0:
            motiflist.append(motif.strip("\""))
    return motiflist

def readTXT():
    '''
    Het bestand met het MSA resultaat
    (langste gen en percentage conservering)
    wordt uitgelezen en geformateerd. Beide waarden worden samen in een
    lijst teruggegeven.
    '''
    openFile = open('MSARESULT.txt', 'r')
    line = openFile.readlines()
    print line
    openFile.close()
    setNames = []
    for x in line:
        entry = x.strip("\n")
        setNames.append(entry.strip("\""))
    return(setNames)

def generatepdf():
    '''
    Een argument parser wordt gemaakt.
    Deze controleert of de benodigde argumenten aanwezig zijn.
    Deze worden opgeslagen in een dictionary
    Alle onderdelen van het verslag worden gegenereerd en meegegeven
    aan een functie die het .tex bestand maakt.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-num','--student_num')
    parser.add_argument('-n','--name')
    parser.add_argument('-t','--title')

    args=parser.parse_args()
    argdict = args.__dict__
    frontpage = front_page(argdict)
    inl = inleiding()
    motifs = readCSV()
    m_en_m = materiaal_methode()
    longGene = readTXT()
    resultList = motifs + longGene
    result = resultaat(resultList)
    discussion_conclusion = discussie_conclusie()
    texfile = frontpage+inl+m_en_m+result+discussion_conclusion
    writefile(texfile)
    
def writefile(texfile):
    '''
    Schrijft het gegenereerde verslag naar een .tex bestand.
    '''
    with open('output.tex','a') as f:
        f.write(texfile)

def main(argv):
    if len(argv) >= 2:
        if argv[1] == "-h" or argv[1] == "--h" or argv[1] == "-help" or argv[1] == "--help":
            showUsageInformation()
        else:
            generatepdf()        
                                 
if __name__ == "__main__":
    main(sys.argv)
    

# Additional information:
# =======================
# Dit script is geschreven door Jesse Kerkvliet en
# kan gebruikt worden bij de opdracht van bapgc door
# Sebastiaan de Vriend.
#
# Het script genereert een .tex file die het verslag kan genereeren.
