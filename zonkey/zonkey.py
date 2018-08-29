from optparse import OptionParser

# parse and extract options from input file
parser = OptionParser()
parser.add_option('-i', '--input', dest='infile', \
                  help='Input file', metavar='FILE')

(options, args) = parser.parse_args()
if not options.infile:   # if filename is not given
    inout.printoptionhelp()
    parser.error('No input file specified\nusage: \'python zonkey.py -i myinputfile\'')

# read input file and obtain a dictionary of options
options = inout.readinputfile(options.infile)


