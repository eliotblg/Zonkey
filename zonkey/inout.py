import os

# read input files and output a dictionary of options

mandatory = ['coordinates']
keywordtype = {'coordinates': 'file', 'qmlist': 'ilist', 'mmlist': 'ilist', \
               'memory': 'float', 'nproc': 'int', 'qmcode': 'string', \
               'qmham': 'string', 'qmbasis': 'string', 'qmextra':'string', \
               'mmcode': 'string', 'mmstruct': 'file', 'mmparam': 'slist', \
               'activelist': 'ilist', 'frozenlist': 'ilist'}
default = {'coordinates': None,'qmlist': [], 'mmlist': [], \
           'memory': 2, 'nproc': 1, 'qmcode': 'gaussian', \
           'qmham': 'hf', 'qmbasis': None, 'qmextra':'', \
           'mmcode': 'namd', 'mmstruct': None, 'mmparam': None, \
           'activelist': [], 'frozenlist': []}
keyworddescription = { \
           'coordinates': 'Coordinate file, can be a .pdb or .xyz', \
           'qmlist': 'List of atoms located in QM region (first atom is 0)', \
           'mmlist': 'List of atoms in MM region, default is reverse of qmlist', \
           }

def readinputfile(ifile):
    finput = default
    if not os.path.exists(ifile) or not os.path.isfile(ifile):
        print('Cannot find input file: ' + ifile)
        exit(1)

    with open(ifile) as f:
        lines = f.readlines()

    for line in lines:
        # read line and take ! and # as comment indicators
        l = line.replace('\=','PLACEHOLDER')
        ll  = l.replace('!','#').split('#')[0] \
               .replace('=',' ').replace('PLACEHOLDER','=').split()
        if len(ll):
            keyword = ll[0]
            if keyword not in keywordtype:
                printoptionhelp()
                print('Unknown keyword: ' + ll[0])
                exit(1)
            if len(ll) < 2:
                printoptionhelp()
                print('Keyword: ' + ll[0] + ' set without value')
                exit(1)
            if keywordtype[keyword][1:] == 'list':
                val = readinputlist(ll)
                if keywordtype[keyword][0] == 'i':
                    val = [int(i) for i in val]
                elif keywordtype[keyword][0] == 'f':
                    val = [float(i) for i in val]
                elif keywordtype[keyword][0] == 's':
                    val = [i for i in val]
            elif keywordtype[keyword] == 'file':
                val = ll[1]
            elif keywordtype[keyword] == 'int':
                val = int(ll[1])
            elif keywordtype[keyword] == 'float':
                val = float(ll[1])
            elif keywordtype[keyword] == 'string':
                val = ll[1].lower()
            finput[keyword] = val    
        
    return finput

# designed to read list directly or read a file from which the list is read
def readinputlist(ilist):
    # if list given directly
    if ilist[1][0] == '[' or ilist[1][0].isdigit():
        llist = ' '.join(ilist[1:]).replace('[','').replace(']','').replace(',',' ').split()
    # if file name is given
    else:
        with open(ilist[1]) as f:
            lines = f.readlines()
        tmplist = [ilist[0]]
        for line in lines:
            tline = line.replace('!','#').split('#')
            if len(tline) and tline[0] != '':
                tmplist.append(tline[0].replace('\n', ''))
        llist = readinputlist(tmplist)
    return llist

def printoptionhelp():
    print('\nYou must provide an input file containing the following keywords:\n')
    print('############################################################')
    print('%15s%15s%15s%15s'%('keyword name','keyword type','default value','mandatory?'))
    print('############################################################')
    for i in keywordtype:
        if i in mandatory:
            print('%15s%15s%15s%15s'%(i,keywordtype[i],default[i],'Yes'))       
        else:
            print('%15s%15s%15s%15s'%(i,keywordtype[i],default[i],'No'))



