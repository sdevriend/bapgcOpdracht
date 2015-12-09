import argparse
import subprocess
import os
import datetime
firstpage=r'''

\documentclass{article}
\usepackage{graphicx}
\begin{document}
\textbf{\Huge %(title)s \\}
\vspace{1cm}
\includegraphics[height=10cm, width=10cm]{Pathwayimage}
\vspace{2cm}
\textbf{\large %(name)s \\}
\vspace{1cm}
\textbf{\normalsize %(student_num)s \\}
\textbf{\normalsize %(month)s \\}
\end{document}
'''
now = datetime.datetime.now()
parser = argparse.ArgumentParser()
parser.add_argument('-num','--student_num')
parser.add_argument('-n','--name')
parser.add_argument('-t','--title')
parser.add_argument("ss")

args=parser.parse_args()
argdict = args.__dict__
month_year = "Pizza"
month_year = "%s %s"%(str(now.strftime("%B")),str(now.year))
print month_year
argdict["month"] = month_year
print args
print args.__dict__
firstpage%argdict

with open('firstpage.tex','w') as f:
    f.write(firstpage%args.__dict__)

os.system('pdflatex firstpage.tex')

os.unlink('firstpage.log')
