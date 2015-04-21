matlab_print_table
==================
Print data in a compact text table format or latex tabular enviroment. 
Can print tables from numeric matrices, cell matrices with different data types or data in the table format.
The second (cell) input specifies the print format (using fprintf syntax), and can be specified for each row, column or for each individual element in the table.
The table is printet with aligned columns and text alignment by inserting spaces. 
Can also print the transposed table version of the supplied data.
Provides a variety of cosmetic options such as adding/removing extra space padding or changing the separation characters.

Syntax:
print_table(dataTable)
print_table(dataCell)

print_table(dataCell, dataDescCellstr)
print_table(dataCell, dataDescCellstr, headerColumnCellstr, headerRowCellstr)

print_table(__, Name, Value)

tableStr = print_table(__)

Input:
dataTable            - table data to print (see table)
dataCell             - cell with data to print (can be numeric matrix)
dataDescCellstr      - cell with sprintf syntax for elements in dataCell       
   Note, dataDesc is expanded to the complete table
   size if only a single element, row or column
   description is supplied. 
headerColumnCellstr  - cell array with column header names
headerRowCellstr     - cell array with row header name

Note, if both headerRow/Column are supplied, one can be one element longer
than the dimension of the dataCell, the extra element (which should be the
first element in the array) is then positioned at the top left corner.

Options, supplied as (..., Name, Value) pairs, overrides default values:
printHeaderCol = 1  - print header columns (if supplied)
printHeaderRow = 1  - print header rows (if supplied)
transposeTable = 0  - transpose table compared to input format
printMode = 'text'  - print mode, 'text' or 'latex'
colSepStr = '|'     - separation string between columns (if 'text')
rowSepStr = ''      - separation line character between rows (if 'text')
rowHSepStr = '-'    - separation line character between header and data
colHSepStr = ''     - extra separation string between col.header and data
textAlignment  = 'c'- text alignment in each column (alt. 'l' or 'r')
  Note, is possible to supply for each column as string, e.g. 'lcl...cr'.
numSpaceColPad = 1  - extra space padding in each column 
spaceColPadAlign = 1- use the extra space padding with the text alignment
  Note, cosmetic change if we do not want the extra space padding to be
  included in the aligned text, e.g. 'lText   ' -> ' lText  ', if false
and numSpaceColPad = 1 and textAlignment = 'l'.
printLatexFull = 1  - add tabular environment to latex table format
printBorder = 0     - print simple border around the table (in text mode)   
borderRowStr = '-'  - border type string, should be single character

Output, table printed in command window, or 
tableStr  - string with output table, preferably printed using fprintf
