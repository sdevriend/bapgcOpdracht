import argparse
import subprocess
import os
import datetime
def front_page(argdict):
    now = datetime.datetime.now()
    firstpage=r'''
    \documentclass{article}
    \usepackage{graphicx}
    \begin{document}
    \textbf{\Huge %(title)s \\}
    \vspace{1cm}
    \begin{figure}[h]
    \centering
    \includegraphics[height=10cm, width=10cm]{test}
    \end{figure}
    \vspace{3cm}
    \textbf{\large %(name)s \\}
    \vspace{1cm}
    \textbf{\normalsize %(student_num)s \\}
    \textbf{\normalsize %(month)s \\}
    \end{document}
    '''
    month_year = "%s %s"%(str(now.strftime("%B")),str(now.year))
    argdict["month"] = month_year
    texfile = firstpage%argdict
    return texfile
def main_function():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-num','--student_num')
    parser.add_argument('-n','--name')
    parser.add_argument('-t','--title')
    parser.add_argument("ss")

    args=parser.parse_args()
    argdict = args.__dict__
    
    if (argdict["ss"] == "FP"):
        texfile = front_page(argdict)
    with open('firstpage.tex','w') as f:
        pass
    writefile(texfile)
    execute()
def writefile(texfile):
    with open('firstpage.tex','a') as f:
        f.write(texfile)

def execute():
    os.system('pdflatex firstpage.tex')


if __name__ == "__main__":
    main_function()

