xls_write_sheets.py
==============================
Automatically save xls files with multiple sheets and lines of 
interest selected (e.g. with significant p-values)


.. automodule:: xls_write_sheets
   :members:
   
   
**List of parameters**

input_type: possible values = 1 or 2. If set to 1, the input is a python dictionary saved in json format with pairs: (names of the Excel sheet, lists of rows to be printed)
e.g: {'sheet1': [header1, row11, row12]: 'sheet2': [header2, row21, row22]}. Headers and rows are lists with strings or numbers; if set to 2, the input is a tab delimited txt files to be converted in a single Excel with multiple sheets. 
In both cases the first line must contain the header, that can be printed in a different format (e.g. bold)

folder: path where the files are

filelist: names of the input files separated by comma with no spaces, or all if entire folder

output_file: name of the resulting xls file

output_folder: name of the results folder

sel_column: index corresponding to the column containing a condition for row selection (e.g. p-values); 0 if none should be selected 

sel_operator: operator for row selection; possible values are > or <

sel_threshold: threshold value for row selection

sel_color: foreground color used for row highlighting



**Examples**

*Python call*

# a) Using the parameters in the parameter file

python xls_write_sheets.py -p xls_parameters_template.txt 
    
# b) Updating some parameters in the parameter file

python xls_write_sheets.py -p xls_parameters_template.txt --input_type 1 --folder myfolder --sel_column 2 --sel_threshold 0.01 --sel_color pink) 

*Calling this script from R*

The R script mesh_xls.R provides an interface with R to both MeSH enrichment and xls writing.
The function xls_write calls the Python script xls_write_sheets.py. An example of use follows:

source('mesh_xls.R')

folder_out <- 'res_xls'

# a) Using the parameters in the parameter file

xls_write(folder_out, par_file='xls_parameters_template.txt')

# b) Updating some parameters in the parameter file

xls_write(folder_out, output='mesh_enriched_all', output_f='test', font='Calibri', sel=2, sel_oper='"<"', sel_thr=0.01, col='pink')

**Notes**

All needed python packages are available on scc, after loading epd/epd-2.7.3 (Traditional Enthought Python Distribution)