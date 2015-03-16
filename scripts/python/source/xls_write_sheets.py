#-------------------------------------------------------------------------------
# Name:        xls_write_sheets
# Purpose:     automatically save xls files with multiple sheets
#              and selected lines (useful for results of enrichment analysis)
#
# Author:      Francesca
#
# Created:     25/03/2014

#-------------------------------------------------------------------------------
#!/usr/bin/env python

"""
Usage: python xls_write.py -p parameter_file.txt
-h help
-p parameter file				*[No default value]
"""

import os
import json
import sys, getopt
import utils

def write_xls_sheets(data, filename, selected=None, fn='Calibri', col_s='yellow'):
    """ writes xls spreadsheeets.
    if selected is provided, rows are highligted accordingly with color = col_s
    :param: data: data instance
    """
    import os
    import xlwt
    from xlwt import Workbook, easyxf
    from xlwt.Utils import rowcol_to_cell
    style_header = easyxf('font: name %s, bold on' %fn)
    stylerow = easyxf('font: name %s' %fn)
    stylerow_selected = easyxf('font: name %s; pattern: pattern solid, fore_colour %s' %(fn, col_s))

    wb = Workbook()
    for dname, dvalues in data.items():
        header = dvalues[0]
        colcount = len(header)
        ws = wb.add_sheet(dname)
        for col in range(colcount):
            ws.write(0, col, str(header[col]), style_header)
            for row in range(1,len(dvalues)):
                if selected and selected[dname][row] == True:
                    s_row = stylerow_selected
                else:
                    s_row = stylerow
                ws.write(row, col, dvalues[row][col], s_row)
    wb.save('%s.xls' %filename)


def check_input(vals):
    """ Takes a list of strings and converts them into floats. """
    newvals = []
    for v in vals:
        try:
            nv = float(v)
        except:
            nv = v
        newvals.append(nv)
    return newvals


def read_files(input_folder, filelist):
    """ Reads the files in filelist from input_folder and saves the corresponding data in a dictionary.
    The function check_input is used to convert string values into floats."""
    dict_files = {}
    for filename in filelist:
        lines = [l.rstrip().split('\t') for l in open(os.path.join(input_folder, filename)).readlines()]
        sheetname = filename.split('.txt')[0]
        dict_files[sheetname] = [check_input(line) for line in lines]

    return dict_files


def main():
    """ Checks parameters and execute write_xls_sheets. If input_type = 1 it reads a json file from
    the folder specified, if input_type = 2 calls read_files """
    list_args = sys.argv[1:]
    print list_args

    if '-h' in list_args:
        print __doc__
        sys.exit(0)
    elif len(list_args) <2 or '-p' not in list_args:
        print __doc__
        sys.exit(0)
    else:
        params, parameter_file, updated = utils.update_parameters(list_args)

        print params

        folder = params['folder']
        filelist_str = params['filelist']
        output_f = params['output_folder']
        if filelist_str == 'all':
           filelist = os.listdir(folder)
        else:
           filelist = filelist_str.split(',')
        if params['input_type'] == '1': # dictionary
            print '%s/%s' %(folder, filelist[0])
            json_data = open('%s/%s' %(folder, filelist[0])).read()
            data = json.loads(json_data)

        elif params['input_type'] == '2': # text files
            data = read_files(folder, filelist)
        else:
            print "error on data type"
            sys.exit(0)

        path = os.getcwd()
        res_folder = os.path.join(path, output_f)
        
        if not(os.path.exists(res_folder)):
            os.makedirs(res_folder)
        outfile = os.path.join(res_folder, params['output_file'])

        if params['sel_column'] != '0':
            selecting_cond = 'x[%s-1]%s%s' %(params['sel_column'], params['sel_operator'], params['sel_threshold'])
            fsel = lambda x: eval(selecting_cond)
            sel = dict([(dname, map(fsel, d)) if len(d[0])>=int(params['sel_column']) else (dname, [False for dd in d]) for dname,d in data.items()])
            write_xls_sheets(data, outfile, sel, fn=params['font'], col_s=params['sel_color'])
        else:
            write_xls_sheets(data, outfile, fn=params['font'])
        print "results saved in %s.xls" %os.path.join(res_folder, params['output_file'])


if __name__ == '__main__':

    main()

