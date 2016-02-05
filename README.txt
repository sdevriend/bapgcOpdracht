Welkom bij de BAPGC-pipeline die met twee simpele opties een pipeline runt
vanaf het begin tot het einde.
De pipeline maakt gebruik van de volgende software packages:
R3.2.2
clustalW2.1
texlive
R heeft de volgende systeem packages nodig:
libcurl4-openssl-dev
libxml2-dev
libgsl0-dev
Binnen R3.2.2 worden de volgende biocLite packages gebruikt:
KEGGREST
GenomicRanges
BSgenome.Hsapiens.UCSC.hg38
org.Hs.eg.db
TxDb.Hsapiens.UCSC.hg38.knownGene
JASPAR2014
TFBSTools
hgu95av2.db
limma
msa
biomaRt
seqinr
ape
ggtree
De packages kunnen gedownload worden door y in te typen bij het startup-script
maar voor het downloaden van de packages moet je wel een root wachtwoord
opgeven omdat het downloadscript met sudo wordt aangeroepen. Als je er
zeker van bent dat de bovenstaande packages aanwezig zijn, dan kan je
n invullen bij de download packages vraag.
Het script kan opgestart worden met:
./StartBapgcPipeline.bash
Hierin wordt gevraagd welke pathway je de analyse voor wilt doen. Voor
de gekozen pathway wordt een resultaat folder en een temp folder gemaakt. De
scripts zijn te vinden in de map Script.

De pipeline en code zijn gemaakt door Sebastiaan de Vriend en Jesse Kerkvliet