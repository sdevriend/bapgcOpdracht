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
    writefile(texfile)
def writefile(texfile):
    with open('output.tex','a') as f:
        f.write(texfile)


if __name__ == "__main__":
    main_function()

