
#-------------------------------------------------------------------------------
# Name:        write doc
# Purpose:     automatically save rst files for sphinx documentation
#
# Author:      Francesca
#
# Created:     25/11/2014

#-------------------------------------------------------------------------------
#!/usr/bin/env python

"""
Usage: python xls_write.py -p parameter_file.txt
-h help
-f script file				*[No default value]
"""

import sys
import os

def add_rst(filename):
    """ add an rst file with docs on each commented function """ 
    file = open('source/%s.rst' %filename, 'w')
    file.write("%s.py\n==============================\n" %filename)
    file.write(".. automodule:: %s\n" %filename)
    file.write("\t:members:")
    file.close()


def update_index(filename):
    """ add the name of the new script to the index file """
    file = open('source/index.rst')

    list_src = []
    flag = 0
    all_lines = []

    for l in file.readlines():
        all_lines.append(l)
        if flag ==1 and set(l) != set([' ', '\n']) and set(l) != set(['\n']) and set(l)!= set(['\r', '\n']) and ("Indices" not in l):
            list_src.append(l)
        if "maxdepth" in l:
            flag = 1
        if "Indices" in l:
            flag = 0

    print list_src
    
    # check if the file is already in the docs, otherwise it doeas not need to be added to the index doc file
    if len([1 for f in list_src if filename in f]) == 0:
       idx = all_lines.index(list_src[-1])
       all_lines.insert(idx+1, "   %s\n" %filename)
       new_index = open("source/index.rst", "w")
       new_index.writelines(all_lines)
       new_index.close()


def main():
    """ Checks parameters and execute the file. """
    list_args = sys.argv[1:]
    print list_args

    if '-h' in list_args:
        print __doc__
        sys.exit(0)
    elif len(list_args) <2 or '-f' not in list_args:
        print __doc__
        sys.exit(0)
    else:
        pnames = list_args[0::2]
        tuple_args = zip(pnames, list_args[1::2])
        file_param = filter(lambda x: x[0]== '-f', tuple_args)
        filen = file_param[0][1]
        file_param_2 = filter(lambda x: x[0]== '-d', tuple_args)
        docn = file_param_2[0][1]
        scripts = os.listdir('source')
        if '%s.rst' %filen in scripts:
           print ".rst file already present in source dir"
        else:
             print "updating rst files"   
             add_rst(filen)
             update_index(filen)
        os.system("python sphinx_mod/sphinx-build.py source %s" %docn)
            
        
if __name__ == '__main__':

    main()

